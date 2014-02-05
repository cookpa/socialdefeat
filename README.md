socialdefeat
============

R package for analysis of social defeat data


## INSTALLATION


## Install dependencies

First, grab the dependencies:

    install.packages(c("cluster", "mclust"), dependencies = TRUE)

If this fails, try setting a CRAN mirror

    chooseCRANmirror()

and then running the installation again.


## Get the source

You can get the source from git or from a zip file.

### Get the source from git

Clone the repository containing the source

    git clone https://github.com/cookpa/socialdefeat.git


### Get the source from a ZIP file

You can grab a ZIP file by going to the project's page and clicking the "Download ZIP" button.

Note that the ZIP file will be called "socialdefeat-master" and will unzip to a folder of that name.
Rename this folder "socialdefeat".



## Install the package

You need to know the path to the directory where you put the source. For example

    sdPath <- "~/code/R/socialdefeat" 

On Windows, you can get the path by holding down shift and right-clicking on the "socialdefeat" folder, 
then selecting "Copy as path".

On a Mac, you can drag the folder from Finder into the R console or a Terminal to see the path.

To check you've got the right path, do

    list.files(sdPath)

You should see

> \> list.files(sdPath)  
> [1] "DESCRIPTION" "man"         "NAMESPACE"   "R"           "README.md"  

Now you are ready to install:

    install.packages(sdPath, repos=NULL, type = "source") 



## Loading the package for use

    library(socialdefeat)


## View the documentation

 
    help(package = socialdefeat)



## Uninstall

First unload the package:

    detach("package:socialdefeat", unload=TRUE)

Then remove it

    remove.packages("socialdefeat")



