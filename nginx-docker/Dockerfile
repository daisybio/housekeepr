FROM nginx

RUN apt-get update
RUN apt-get install -y software-properties-common
RUN apt-get install -y curl
RUN apt-get install -y gpg
RUN curl http://nginx.org/keys/nginx_signing.key | apt-key add -
RUN echo "deb http://nginx.org/packages/debian/ stretch nginx\n deb-src http://nginx.org/packages/debian/ stretch nginx" >> /etc/apt/sources.list
RUN cat /etc/apt/sources.list
RUN apt-get update
RUN apt-get install -y dpkg-dev
WORKDIR /tmp
RUN apt-get source nginx
RUN apt-get build-dep -y nginx
RUN mv nginx-* nginx
WORKDIR /tmp/nginx
RUN apt-get install -y git
RUN git clone https://github.com/thomsonreuters/nginx-sticky-module-ng
RUN rm nginx-sticky-module-ng/docs/sticky.pdf
RUN rm nginx-sticky-module-ng/docs/sticky.vsd

RUN LINE=`grep -n "config.status.nginx:" debian/rules | cut -f1 -d:` && LINE=$(($LINE+2)) && sed -i $LINE's/$/ --add-module=\/tmp\/nginx\/nginx-sticky-module-ng/' debian/rules
RUN sed -i 's/$CFLAGS \-Werror//' ./auto/cc/icc ./auto/cc/gcc ./auto/cc/clang 
RUN EDITOR=/bin/true dpkg-source -y --commit ./ sticky
RUN dpkg-buildpackage -uc -us

WORKDIR /tmp
RUN dpkg --install nginx_*.deb 

COPY nginx.conf /etc/nginx/nginx.conf
ADD ./wait-for-it.sh  /
ADD ./www /www/

