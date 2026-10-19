# syntax=docker/dockerfile:1.7
FROM centos:6

COPY CentOS-Vault.repo /etc/yum.repos.d/CentOS-Vault.repo
RUN find /etc/yum.repos.d -name 'CentOS-*.repo' ! -name 'CentOS-Vault.repo' -delete \
    && yum install -y tcl tk xterm perl which libX11 libstdc++ \
    && yum clean all

# The bind mount keeps the IBM RPMs out of Docker's build layers. The image
# contains only their installed payload, and is intended for local use.
RUN --mount=type=bind,from=cell_rpms,source=systemsim-cell-3.1-25.f9.x86_64.rpm,target=/tmp/systemsim.rpm,ro \
    --mount=type=bind,from=cell_rpms,source=sysroot_image-3.1-1.noarch.rpm,target=/tmp/sysroot.rpm,ro \
    rpm -ivh --nodeps --noscripts /tmp/systemsim.rpm /tmp/sysroot.rpm \
    && test -x /opt/ibm/systemsim-cell/bin/systemsim-cell \
    && test -s /opt/ibm/systemsim-cell/images/cell/sysroot_disk \
    && test -s /opt/ibm/systemsim-cell/images/cell/vmlinux

ENV PATH=/opt/ibm/systemsim-cell/bin:$PATH
RUN mkdir -m 1777 /sim
WORKDIR /sim
CMD ["systemsim", "-cell"]
