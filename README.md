# Gemini Reductions for ZTF

Use the scripts here to reduce Gemini spectra for ZTF infant SNe set-up

## Software required

These instructions will get you a copy of the project up and running on your local machine. Install Astroconda python virtual environment following the instructions here:
https://astroconda.readthedocs.io/en/latest/getting_started.html

Make sure to get the version with iraf and python 2.7.X

On Mac OS X you also need XQuartz.


Grab the task LA Cosmic for cosmic ray rejection, instructions for doing so can be found on Gemini's own data reduction cookbook:
http://ast.noao.edu/sites/default/files/GMOS_Cookbook/GettingStarted.html
or try this link:
http://www.astro.yale.edu/dokkum/lacosmic/lacos_spec.cl

### Set-up

After installation and set up start your virtual environment (mine is called iraf27) then issue the command mkiraf in the directory you will be working in.

```
source activate iraf27
mkiraf
```

Afterwards edit the login.cl file to contain the lacos_spec.cl file as a task using the following command

```
task lacos_spec= /path/to/iraf_home/lacos_spec.cl
```

### Data & reduction

Grab the data from archive.gemini.edu

Make sure to unzip you data such as by using

```
gunzip *.bz2
```

Then we run the script to organize the data:

```
python readsort.py
```
Check to make sure you have all the data at this point. A bunch of text files linking to the fits files should be created with names such as object_b600_520.txt
Finally run the reduction script.

```
python ztf_reduce.py
```

Which will guide you through a semi-automated standard IRAF reduction. You will have to check and confirm the steps as they are shown. For more details or finer control, see the cookbook with the step by step instructions [Gemini Cookbook V2.docx](Gemini Cookbook V2.docx)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Disclaimer: This is barely used minimally functioning code.
I basically used it to get rapid first look reductions for iPTF (and now ZTF) Gemini spectra.

This code was written *years ago* and sparsely updated to simply function at a bare minimum level. The original goal was to take a series of IRAF reduction steps and transparently put them into a customizable python script. If people end up using it let me know and I can generalize it. If you notice something wrong please also email me. You can [find my contact information here](http://www.su.se/english/profiles/emka6994). Absolutely no one should use this to judge the level of my coding skill or style.

## Acknowledgements
Everything here is based on a workshop given by Brad Cenko at an iPTF summer school in 2014. Christoffer Fremling & Emir Karamehmetoglu went through the initial drafting phase in 2014/2015. Since then it has been maintained and updated by Emir Karamehmetoglu several times.
