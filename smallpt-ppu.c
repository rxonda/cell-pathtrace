#define _XOPEN_SOURCE 600
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <libspe2.h>
#include <pthread.h>
#include <altivec.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <string.h>
#define MAX_SPU (6)
#ifndef WIDTH
#define WIDTH (640)
#endif
#ifndef HEIGHT
#define HEIGHT (480)
#endif


//char * spu_program = "./smallptSPU";
//char * spu_program = "./smallptSPUdbl";
//char * spu_program = "./smallptSPUdbl2";
//char * spu_program = "./smallptSPUdbl3";
//char * spu_program = "./forwardSPU";
char * spu_program = "./forwardFloatSPU";

typedef vector signed int Vecint;
Vecint vecint(int x, int y, int z){
	Vecint r={x,y,z,0};
	return r;
}

typedef struct _contexto {
	unsigned long long retorno, status_termination, status; //endereço de retorno
	unsigned long long seed;
	int samps, w, h, dummy;
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
int worker_index[MAX_SPU];

int spus=MAX_SPU, samps=1;
unsigned long long render_seed=1;

int pathtrace(int segment) {
	spe_stop_info_t stop={0};
	unsigned int entry = SPE_DEFAULT_ENTRY;
	int ret;

	ctxt[segment].w=WIDTH;
	ctxt[segment].h=HEIGHT;
	ctxt[segment].samps = samps;
	ctxt[segment].seed = render_seed;
	ctxt[segment].status_termination=(unsigned long long)&status_termination[segment];
	ctxt[segment].status=(unsigned long long)status;
	ctxt[segment].retorno=(unsigned long long)c_int;

	ret=spe_context_run(spe[segment], &entry, 0, &ctxt[segment], NULL, &stop);

	if(ret<0 || stop.stop_reason!=SPE_EXIT || stop.result.spe_exit_code!=0) {
		fprintf(stderr,"SPE %d failed (run=%d, reason=%u, exit=%d)\n",
		        segment,ret,stop.stop_reason,stop.result.spe_exit_code);
		return -1;
	}
	return 0;
}

void *local_thread(void *data) {
	int segment=*(int*)data;
	return pathtrace(segment)==0 ? NULL : (void*)1;
}

int run_spu(void) {
	spe_program_handle_t * prog;
	int ret, i;

	prog=spe_image_open(spu_program);

	if(!prog) {
		fprintf(stderr,"Cannot open SPU image: %s\n",spu_program);
		return -1;
	}

	for(i=0; i<spus; i++) {
		spe[i]=spe_context_create(0, NULL);
		if(!spe[i]) {
			fprintf(stderr,"Failed to create SPE context %d\n",i);
			return -1;
		}
		ret=spe_program_load(spe[i], prog);
		if(ret) {
			fprintf(stderr,"Failed to load SPU program on SPE %d\n",i);
			return -1;
		}
		worker_index[i]=i;
		ret=pthread_create(&spuid[i], NULL, local_thread, &worker_index[i]);
		if(ret) {
			fprintf(stderr,"Failed to start SPE thread %d: %s\n",i,strerror(ret));
			return -1;
		}
	}
	spe_image_close(prog);
	return 0;
}

int main(int argc, char** argv) {
	int w=WIDTH, h=HEIGHT, i;
	FILE *f;
	char *end;
	long value;
	int failed=0;
	if (argc > 4) {
		fprintf(stderr,"Usage: %s [samples: multiple of 4] [SPEs: 1..6] [seed]\n",argv[0]);
		return 2;
	}
	if (argc >= 2) {
		errno=0;
		value=strtol(argv[1],&end,10);
		if(errno || *end || value<4 || value%4 || value>INT_MAX) {
			fprintf(stderr,"Samples must be a positive multiple of four\n");
			return 2;
		}
		samps=(int)value/4;
	}
	if (argc >= 3) {
		errno=0;
		value=strtol(argv[2],&end,10);
		if(errno || *end || value<1 || value>MAX_SPU) {
			fprintf(stderr,"SPE count must be between 1 and %d\n",MAX_SPU);
			return 2;
		}
		spus=(int)value;
	}
	if (argc >= 4) {
		errno=0;
		render_seed=strtoull(argv[3],&end,0);
		if(errno || end==argv[3] || *end || argv[3][0]=='-') {
			fprintf(stderr,"Seed must be an unsigned 64-bit integer\n");
			return 2;
		}
	}
	printf("Starting %d SPE worker(s), %d samples per pixel, seed %llu\n",
	       spus,samps*4,render_seed);
	status[0]=0;
	if (run_spu()) return 1;
	for (i=0; i<spus; i++) {
		void *result=NULL;
		int ret=pthread_join(spuid[i],&result);
		if(ret || result) {
			fprintf(stderr,"SPE worker %d failed\n",i);
			failed=1;
		}
		if(spe_context_destroy(spe[i])) {
			fprintf(stderr,"Failed to destroy SPE context %d\n",i);
			failed=1;
		}
	}
	if(failed) return 1;

	f=fopen("image.ppm","w");
	if(!f) {
		perror("image.ppm");
		return 1;
	}
	if(fprintf(f,"P3\n%d %d\n255\n",w,h)<0) failed=1;
	for(i=0; i<w*h && !failed; i++) {
		int pixel[4];
		memcpy(pixel,&c_int[i],sizeof(pixel));
		if(fprintf(f,"%d %d %d ",pixel[0],pixel[1],pixel[2])<0)
			failed=1;
	}
	if(fclose(f)) failed=1;
	if(failed) {
		perror("Failed to write image.ppm");
		return 1;
	}
	printf("\nRender complete\n");
	return 0;
}
