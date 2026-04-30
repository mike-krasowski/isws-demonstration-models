FROM nginx:1.10.1-alpine
COPY src/html /usr/share/nginx/html

#documentation
EXPOSE 80:80

#CMD ["ngnix" , "-g", daemon off ;"]
