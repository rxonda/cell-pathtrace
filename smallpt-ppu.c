#define _XOPEN_SOURCE 600
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <libspe2.h>
#include <pthread.h>
#include <altivec.h>
#define MAX_SPU (6)
#define WIDTH (640)
#define HEIGHT (480)
#define BUFFERSIZE (1<<14)


//char * spu_program = "./smallptSPU";
//char * spu_program = "./smallptSPUdbl";
//char * spu_program = "./smallptSPUdbl2";
//char * spu_program = "./smallptSPUdbl3";
//char * spu_program = "./forwardSPU";
char * spu_program = "./forwardFloatSPU";

typedef vector int Vecint;
Vecint vecint(int x, int y, int z){
	Vecint r={x,y,z,0};
	return r;
}

typedef struct _contexto {
	unsigned long long retorno, status_termination, status; //endereço de retorno
	int samps, w, h, dummy[3]; //32bytes
} contexto;

typedef struct {
	int status __attribute__((aligned(16)));
} spe_remote_function_status_t;

contexto ctxt[MAX_SPU] __attribute__((aligned(16)));

Vecint c_int[WIDTH*HEIGHT] __attribute__((aligned(128)));

volatile unsigned int status[128/sizeof(unsigned int)] __attribute__((aligned(128)));

spe_remote_function_status_t status_termination[MAX_SPU] __attribute__((aligned(16)));

spe_context_ptr_t spe[MAX_SPU];

pthread_t spuid[MAX_SPU];

int spus=MAX_SPU, samps=1;

void pathtrace(int segment) {
	spe_stop_info_t stop;
	unsigned int entry = SPE_DEFAULT_ENTRY;
	int ret;

	ctxt[segment].w=WIDTH;
	ctxt[segment].h=HEIGHT;
	ctxt[segment].samps = samps;
	ctxt[segment].status_termination=(unsigned long long)&status_termination[segment];
	ctxt[segment].status=(unsigned long long)status;
	ctxt[segment].retorno=(unsigned long long)c_int;

	ret=spe_context_run(spe[segment], &entry, 0, &ctxt[segment], NULL, &stop);

	if(ret<0) {
		printf("Erro ao executar o programa na SPE!\n");
		exit(1);
	}
}

void local_thread(void *data) {
	pathtrace((int)data);
}

typedef struct _Data_t {
	unsigned long long ptr;
	struct _Data_t* anterior;
} Data_t;

typedef struct _Pilha_t {
	Data_t *topo;
	int size;
} Pilha_t;

void push(Pilha_t* pilha,unsigned long long value) {
	Data_t *temp=(Data_t*)malloc(sizeof(Data_t));
	temp->ptr=value;
	temp->anterior=pilha->topo;
	pilha->size++;
	pilha->topo=temp;
}

Data_t pop(Pilha_t* pilha) {
	Data_t retorno, *temp;
	if(pilha->size>0) {
		pilha->size--;
		retorno=*pilha->topo;
		temp=pilha->topo;
		pilha->topo=pilha->topo->anterior;
		free(temp);
	}
	return retorno;
}

Pilha_t pilha[MAX_SPU];

void thread_pilha(void* data) {
	int i=(int)data;
	do {
		unsigned int msg;
		unsigned int status;
		spe_out_intr_mbox_read(spe[i], &msg, 1, SPE_MBOX_ALL_BLOCKING);
		if(msg==65) { //signaling data retrieving
			unsigned long long ptr;
			Data_t temp;
			temp=pop(&pilha[i]);
			ptr=temp.ptr;
			unsigned int ea_hi=(unsigned int)(ptr>>32);
			spe_in_mbox_write(spe[i],&ea_hi, 1, SPE_MBOX_ALL_BLOCKING);
			unsigned int ea_lo=(unsigned int)ptr;
			spe_in_mbox_write(spe[i],&ea_lo, 1, SPE_MBOX_ALL_BLOCKING);
			spe_out_intr_mbox_read(spe[i], &msg, 1, SPE_MBOX_ALL_BLOCKING);
			if(msg==75)
				free(ptr); //freeing some space, wainting for the transfer to complete
		}
		if(msg==55) {
			unsigned long long ptr;
			if(posix_memalign((void*)&ptr,128,BUFFERSIZE)!=0)
				printf("Falhou a alocação de memória da pilha\n");
			push(&pilha[i],ptr);
			unsigned int ea_hi=(unsigned int)(ptr>>32);
			spe_in_mbox_write(spe[i],&ea_hi, 1, SPE_MBOX_ALL_BLOCKING);
			unsigned int ea_lo=(unsigned int)ptr;
			spe_in_mbox_write(spe[i],&ea_lo, 1, SPE_MBOX_ALL_BLOCKING);
		}
	} while(1);
}

void wait_spu(void) {
	int i;
	for(i=0; i<spus; i++)
		while(1)
			if(status_termination[i].status==1)
				break;
}

void run_spu() {
	spe_program_handle_t * prog;
	int ret, i;

	prog=spe_image_open(spu_program);

	if(!prog) {
		printf("Impossivel abrir imagem!\n");
		exit(1);
	}

	for(i=0; i<spus; i++) {
		spe[i]=spe_context_create(0, NULL);
		if(!spe[i]) {
			printf("Erro ao criar contexto da SPE!\n");
			exit(1);
		}
		ret=spe_program_load(spe[i], prog);
		if(ret) {
			printf("Erro ao carregar programa!\n");
			exit(1);
		}
		pthread_create(&spuid[i], NULL, local_thread, (void*)i);
		pthread_create(&spuid[i], NULL, thread_pilha, (void*)i); //Thread controller of the SPE's stack
	}
}

int main(int argc, char** argv) {
	char * filename;
	/*vector unsigned int ccc={10,10,10,10};
	vector unsigned int temp;
	unsigned int * eacc=&temp;
	temp=vec_splats((unsigned int)vec_extract((vector unsigned char)temp,0x4));
	fprintf(stderr, "Teste: %d %d %d %d\n",eacc[0],eacc[1],eacc[2],eacc[3]);
	return 0;*/
	int w=WIDTH, h=HEIGHT;
	if(argc>=2)
		samps=atoi(argv[1])/4;

	if(argc>=3){
		spus=atoi(argv[2]);
		if(spus>MAX_SPU) spus=MAX_SPU;
	}
	printf("Chamando as spus\n");
	status[0]=0; //linha zero
	run_spu();

	wait_spu(); //espera o calculo de cada spu

	FILE *f = fopen("image.ppm", "w");         // Write image to PPM file.
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
	int i;
	for (i=0; i<w*h; i++){
		int *tmp=&c_int[i];
		fprintf(f,"%d %d %d ", tmp[0], tmp[1], tmp[2]);
		
	}
	printf("\nExecucao terminada com sucesso!\n");

	return 0;
}
