import sys
import os
import astropy.io.fits as fits
import glob

#logincl_dir='/Users/emir/iraf'
#workingdir=os.getcwd()
#os.chdir(logincl_dir)
from pyraf import iraf
#os.chdir(workingdir)


iraf.gemini()
iraf.gmos()
iraf.stsdas()
iraf.onedspec()

yes='yes'
no='no'
default_caldir='onedstds$iidscal/'


l=open('reduction.log','w')

#Produce normalized flats
#if not os.path.isfile('b600_525_norm_flat_3.fits'):
iraf.gsflat(inflats="@flatb600_520_1.txt", specflat='b600_520_norm_flat_1.fits', order=20, fl_over=yes, fl_bias=no, fl_inter=no, fl_detec=yes, fl_seprows=no)
iraf.gsflat(inflats='@flatb600_520_2.txt', specflat='b600_520_norm_flat_2.fits', order=20, fl_over=yes, fl_bias=no, fl_inter=no, fl_detec=yes, fl_seprows=no)
iraf.gsflat(inflats='@flatb600_520_3.txt', specflat='b600_520_norm_flat_3.fits', order=20, fl_over=yes, fl_bias=no, fl_inter=no, fl_detec=yes, fl_seprows=no)
iraf.gsflat(inflats='@flatb600_525_1.txt', specflat='b600_525_norm_flat_1.fits', order=20, fl_over=yes, fl_bias=no, fl_inter=no, fl_detec=yes, fl_seprows=no)
iraf.gsflat(inflats='@flatb600_525_2.txt', specflat='b600_525_norm_flat_2.fits', order=20, fl_over=yes, fl_bias=no, fl_inter=no, fl_detec=yes, fl_seprows=no)
iraf.gsflat(inflats='@flatb600_525_3.txt', specflat='b600_525_norm_flat_3.fits', order=20, fl_over=yes, fl_bias=no, fl_inter=no, fl_detec=yes, fl_seprows=no)

#log?
l.write('flats=yes')

#Reduce Science Arc and Standard Star Images


iraf.gsreduce(inimages='@object_b600_520_1.txt', flatim='b600_520_norm_flat_1.fits', fl_bias=no, fl_crspec=yes)
iraf.gsreduce(inimages='@object_b600_520_2.txt', flatim='b600_520_norm_flat_2.fits', fl_bias=no, fl_crspec=yes)
iraf.gsreduce(inimages='@object_b600_520_3.txt', flatim='b600_520_norm_flat_3.fits', fl_bias=no, fl_crspec=no)

iraf.gsreduce(inimages='@object_b600_525_1.txt', flatim='b600_525_norm_flat_1.fits', fl_bias=no, fl_crspec=yes)
iraf.gsreduce(inimages='@object_b600_525_2.txt', flatim='b600_525_norm_flat_2.fits', fl_bias=no, fl_crspec=yes)
iraf.gsreduce(inimages='@object_b600_525_3.txt', flatim='b600_525_norm_flat_3.fits', fl_bias=no, fl_crspec=no)

#if not os.path.isfile('b600_520_norm_flat_1.fits'):

iraf.gsreduce(inimages='@CuAr_b600_520_1.txt', flatim='b600_520_norm_flat_1.fits', fl_bias=no, fl_crspec=no)
iraf.gsreduce(inimages='@CuAr_b600_525_1.txt', flatim='b600_525_norm_flat_1.fits', fl_bias=no, fl_crspec=no)


#log?
l.write('reduced=yes')

#Wavelength calibrate ARCs in 2d

iraf.gswavelength(inimages='gs//@CuAr_b600_520_1.txt')
iraf.gswavelength(inimages='gs//@CuAr_b600_525_1.txt')

#log?
l.write('arcs=yes')

#Applying Wavelength Calibration to Science Images

iraf.gstransform(inimages='gs//@object_b600_520_1.txt', wavtran='gs//@CuAr_b600_520_1.txt')
iraf.gstransform(inimages='gs//@object_b600_520_2.txt', wavtran='gs//@CuAr_b600_520_1.txt')
iraf.gstransform(inimages='gs//@object_b600_520_3.txt', wavtran='gs//@CuAr_b600_520_1.txt')

iraf.gstransform(inimages='gs//@object_b600_525_1.txt', wavtran='gs//@CuAr_b600_525_1.txt')
iraf.gstransform(inimages='gs//@object_b600_525_2.txt', wavtran='gs//@CuAr_b600_525_1.txt')
iraf.gstransform(inimages='gs//@object_b600_525_3.txt', wavtran='gs//@CuAr_b600_525_1.txt')

#log?
l.write('wavelength=yes')


#Function to read SN
def read_sn(filename):
	f=open(filename,'r')
	obje=f.readline()[:-1]
	f.close()
	return obje

#Function to check if SN or CompStar
def issn(obje):
	import astropy.io.fits as fits
	hdulist= fits.open(obje)
	hdulist.close()
	if 'ZTF' in str(hdulist[0].header['object']):
		print str(hdulist[0].header['object'])+' is a SN'
		return 'ztfsn'
	else:
		print str(hdulist[0].header['object'])+' is a Comp Star'
		return 'ztfstar'

#Function to read comp star data in IRAF format
def read_comp(fname):
	import astropy.io.fits as fitse
	hdulist=fits.open(fname)
	hdulist.close()
	#name in iraf format
	cname=str(hdulist[0].header['object']).translate(None,'+-_').lower()
	amass=float(str(hdulist[0].header['airmass']))
	expotime=float(str(hdulist[0].header['exptime']))
	return cname,amass,expotime

objectlist=glob.glob('object_b600*.txt')

#put all fits files into this dictionary
objects_dict=dict()
objects_pars=dict()
imcopy_names=dict()
for obj_file in objectlist:
	objects_dict[obj_file]=read_sn(obj_file)

	ending=obj_file.split('_')
	tempend='_'+ending[-2]+'_'+ending[-1].split('.')[0]

	tempout=issn(objects_dict[obj_file])+tempend


	ptfsn1_name,ptfsn1_air,ptfsn1_exp=read_comp(objects_dict[obj_file])

	if 'sn' in tempout:
		objtype='sn'
	else:
		objtype='star'
	objects_pars[tempout]=(ptfsn1_name,ptfsn1_air,ptfsn1_exp,objtype)
	obj_TGS_SCI='tgs'+objects_dict[obj_file]+'[SCI]'

	iraf.imcopy(input=obj_TGS_SCI,output=tempout)

	imcopy_names[obj_file]=tempout

#Reduce with apall
for imname in objects_pars.keys():
	print(imname, imname+'ms')
	templlist=list(objects_pars[imname])
	templlist.append(imname+'ms')
	objects_pars[imname]=templlist


for imname in objects_pars.keys():
	iraf.apall(input=imname, out=imname+'ms', nfind=1, interactive=yes, find=yes, recenter=yes, resize=yes, edit=yes, trace=yes, fittrace=yes, t_nsum = 50, t_nlost = 50, extract=yes, extras=yes, review=yes, background='fit')

l.write('extracted=yes')

for imname in objects_pars.keys():
	if 'ztfstar' in imname:
		print(objects_pars[imname][4])
		iraf.standard(input=objects_pars[imname][4], output='std_'+imname, star_name=objects_pars[imname][0], airmass=objects_pars[imname][1], exptime=objects_pars[imname][2], answer=yes,caldir=default_caldir, observatory='Keck')

for imname in objects_pars.keys():
	if 'ztfstar' in imname:
		iraf.sensfunc(standards='std_'+imname, sensitivity='sens_'+imname, observatory='Keck', order=4)
		if '520' in imname:
			sensstar520='sens_'+imname
		elif '525' in imname:
			sensstar525='sens_'+imname


combname=[]
for imname in objects_pars:
	if 'sn' in objects_pars[imname][3]:
		print('name,airmass,exptime,obstype')
		print(objects_pars[imname])
		if '520' in imname:
			sensstar520='sens_ztfstar_520_3.fits'
			usesens=sensstar520
		elif '525' in imname:
			sensstar525='sens_ztfstar_520_3.fits'
			usesens=sensstar525
		iraf.onedspec.calibrate(input=objects_pars[imname][4], out=imname+'_calib', extinction='onedstds$kpnoextinct.dat', observatory='Keck', ignoreaps=yes, sensitivity=usesens,airmass=objects_pars[imname][1],exptime=objects_pars[imname][2])
		combname.append(imname+'_calib')

combstr=','.join(combname)

combname=[]
for imname in objects_pars:
	if 'sn' in objects_pars[imname][3]:
		combname.append(imname+'_calib')
combstr=','.join(combname)
iraf.scombine(input=combstr, out='combinedspec_rmlow.ms', combine='average',reject='minmax', scale='none', sample='5000:7500',nlow=1,nhigh=1,nkeep=2)
iraf.scombine(input=combstr, out='combinedspec.ms', combine='average',reject='minmax', scale='none', sample='5000:7500',nlow=0,nhigh=1,nkeep=1)


l.write('calib=yes')

iraf.splot(images='combinedspec_rmlow.ms.fits')

for objname in objects_pars:
	if 'sn' in objname:
		SNName=objects_pars[objname][0]
		break


iraf.wspectext(input='combinedspec.ms.fits[*,1]',output=SNName+'.ascii',header='NO')
iraf.wspectext(input='combinedspec_rmlow.ms.fits[*,1]',output=SNName+'_rmlow.ascii',header='NO')

l.write('finalname='+SNName+'.ascii')
l.close()
