FROM julia:1.4.1
COPY . /OhMyQSIM
RUN julia -e "using Pkg; Pkg.add(\"OhMyQSIM\")" 
