
This file provide an overview of all the files that you can find in this
project. Also, in this document you can find how the folder structure is set up
and how to start running the program.

THE FOLDER ORGANIZATION GOES AS FOLLOWS: 

	/doc -- Documentation of the project, graphics, etc....

	/------ /verbose.txt -- A lengthy explanation on how the main R script
                                works and how does it connect to C++. How the
                                results are stored. Design decisions, etc...

	/res -- Resources for the project.

	/------ /Cas9_toy    -- An example of the reads. Contain a lightweight
                                real dataset. Excellent to test that your
                                program works in less than a minute.

	/------ /HTML        -- Contain common files that are use in the final
                                HTML report. Such as the .CSS file, the
                                JavaScript file and the logo for the top
                                header.

	/------ /Results_toy -- This is an example of what would happen if you
                                run the content of Cas9_toy with the default
                                parameters.

	/Rpackages -- The packages generated for CRAN. You will find duplicates
                      of the source files here. The folders are set up so you
                      just need to build and check the package in R. You
                      shouldn't get any errors, warnings, or notes.

	/src -- All kind of source codes.

	/------ /CPP    -- The main c++ code, object, makefiles, binaries,
                           etc... goes in here. In particular this project
                           implement the Gotoh sequene alignment algorithm in
                           C++. You can find all of that in here.

	/------ /drafts -- Files that are messy; nice place to store daily
                           backups or where to summit unstable code.

	/------ /Python -- We have a program that correct errors in Fastq files
                           regarding how they are organized with respect each
                           barcode. This problem is call $NAME and you can find
                           the appropiate source code to run that program here.

	/------ /R      -- The main R scripts go here. You will find 4 main
                           scripts. However, these scripts uses more libraries
                           which contain many functions. Those can be found in
                           the library folder.

	/------ /----- /libraries -- The functions used in the main R script
                                     are located here.

	/------ /----- /Rpackages -- We build these packages based on all the
                                     source code mentiened before. In here we
                                     include a very important package called
                                     "gotoh". This is impportant because is the
                                     algorithm to align sequences, and needs to
                                     be into a package because of how the
                                     parallel library "doParallel" works.

	/test -- Tests to check the validity of the code after an update.

	/tutorials -- Tiny pieces of code that show you how to work with
                      certain libraries. In here we also include benchmarks for
                      the different parts and to compare with other similar
                      libraries out there.

	README.txt This document.
