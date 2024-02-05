# Use the official Julia image as the base image
FROM julia:latest
# FROM consumere/shinyapps:japp

# Set the working directory inside the container
WORKDIR /app

# Copy the Julia Project.toml and Manifest.toml files to the container
# COPY Project.toml Manifest.toml /app/ #<-The Manifest.toml is a snapshot of the exact state of a Julia environment. I

COPY Project.toml /app/

#RUN julia -e 'using Pkg; Pkg.add.(["CSV", "DataFrames", "Dash","Plots", "PlotlyJS","Base64","Dates","Statistics"])'
#RUN julia -e 'using Pkg; Pkg.add.(["CSV", "DataFrames", "Dash", "PlotlyJS","Base64","Dates","Statistics"])'
#RUN julia -e 'using Pkg; Pkg.add.(["WaSiM"])' #its dev, so... nope.

# Install the Julia packages
RUN julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate();'
# RUN julia -e 'using Pkg; Pkg.activate("."); Pkg.update(); Pkg.instantiate();'
# Copy the rest of the application files to the container
COPY . /app

#EXPOSE 8080

#/usr/local/bin/docker-entrypoint.sh: 11: exec: julia --optimize=3 --math-mode=fast: not found 
#CMD ["julia", "--math-mode=fast","dpr_dash.jl"] #das funkt nicht wg dem script.
# CMD ["julia", "appts.jl"]
# CMD ["julia", "dpr_dash.jl"]

#CMD ["julia","--math-mode=fast","--optimize=3", "appts.jl"]
#"--optimize=3": The --optimize flag allows you to set the level of optimization. 
#Level 3 is the highest optimization level, which should help reduce memory consumption. 
#It optimizes the code aggressively but might increase compilation time.

#docker build -t consumere/shinyapp:dprjs .
#docker run -p 8080:8050 consumere/shinyapp:dprjs
#docker build -t consumere/shinyapps:dprjs .
#docker run -p 8050:8050 consumere/shinyapps:dprjs