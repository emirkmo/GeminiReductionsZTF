import sys 
from pyraf import iraf 

iraf.gemini()
iraf.gmos()
iraf.stsdas()

yes='yes'
no='no'
default_caldir='onedstds$iidscal/'


l=open('reduction.log','w')

#Produce normalized flats

iraf.gsflat(inflats='@flatb600_1.txt',specflat='b600_norm_flat_1.fits',fl_bias=no,fl_inter=no,fl_detec=yes,fl_seprows=no)

iraf.gsflat(inflats='@flatr400_1.txt',specflat='r400_norm_flat_1.fits',fl_bias=no,fl_inter=no, fl_detec=yes, fl_seprows=no)

iraf.gsflat(inflats='@flatb600_2.txt', specflat='b600_norm_flat_2.fits', fl_bias=no, fl_inter=no, fl_detec=yes, fl_seprows=no)

iraf.gsflat(inflats='@flatr400_2.txt', specflat='r400_norm_flat_2.fits', fl_bias=no, fl_inter=no, fl_detec=yes, fl_seprows=no)

#log?
l.write('flats=yes')

#Reduce Science Arc and Standard Star Images

iraf.gsreduce(inimages='@object_r400_1.txt', flatim='r400_norm_flat_1.fits', fl_bias=no)
iraf.gsreduce(inimages='@CuAr_r400_1.txt', flatim='r400_norm_flat_1.fits', fl_bias=no)
iraf.gsreduce(inimages='@object_b600_1.txt', flatim='b600_norm_flat_1.fits', fl_bias=no)
iraf.gsreduce(inimages='@CuAr_b600_1.txt', flatim='b600_norm_flat_1.fits', fl_bias=no)
iraf.gsreduce(inimages='@object_r400_2.txt', flatim='r400_norm_flat_2.fits', fl_bias=no)
iraf.gsreduce(inimages='@object_b600_2.txt', flatim='b600_norm_flat_2.fits', fl_bias=no)
iraf.gsreduce(inimages='@CuAr_r400_2.txt', flatim='r400_norm_flat_2.fits', fl_bias=no) 
iraf.gsreduce(inimages='@CuAr_b600_2.txt', flatim='b600_norm_flat_2.fits', fl_bias=no)

#log?
l.write('reduced=yes')

#Wavelength calibrate ARCs in 2d 
iraf.gswavelength(inimages='gs//@CuAr_r400_1.txt')
iraf.gswavelength(inimages='gs//@CuAr_b600_1.txt')
iraf.gswavelength(inimages='gs//@CuAr_r400_2.txt')
iraf.gswavelength(inimages='gs//@CuAr_b600_2.txt')

#log?
l.write('arcs=yes')

#Applying Wavelength Calibration to Science Images
iraf.gstransform(inimages='gs//@object_r400_1.txt', wavtran='gs//@CuAr_r400_1.txt')
iraf.gstransform(inimages='gs//@object_b600_1.txt', wavtran='gs//@CuAr_b600_1.txt')
iraf.gstransform(inimages='gs//@object_b600_2.txt', wavtran='gs//@CuAr_b600_2.txt')
iraf.gstransform(inimages='gs//@object_r400_2.txt', wavtran='gs//@CuAr_r400_2.txt')

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
	if 'PTF' in str(hdulist[0].header['object']):
		print str(hdulist[0].header['object'])+' is a SN'
		return 'ptfsn'
	else:
		print str(hdulist[0].header['object'])+' is a Comp Star'
		return 'ptfstar'
		
#Function to read comp star data in IRAF format
def read_comp(fname):
	import astropy.io.fits as fits
	hdulist=fits.open(fname)
	hdulist.close()
	#name in iraf format
	cname=str(hdulist[0].header['object']).translate(None,'+-_').lower()
	amass=float(str(hdulist[0].header['airmass']))
	expotime=float(str(hdulist[0].header['exptime']))
	return cname,amass,expotime	

#Read the Gemini outputs
obj_r1=read_sn('object_r400_1.txt')
obj_r2=read_sn('object_r400_2.txt')
obj_b1=read_sn('object_b600_1.txt')
obj_b2=read_sn('object_b600_2.txt')

#Redparts
out_r1=issn(obj_r1)+'1.fits'

if 'sn' in out_r1:
	ptfsn1_name,ptfsn1_air,ptfsn1_exp=read_comp(obj_r1)
else:
	ptfstar1_name,ptfstar1_air,ptfstar1_exp=read_comp(obj_r1)

out_r2=issn(obj_r2)+'1.fits'

if 'sn' in out_r2:
	ptfsn1_name,ptfsn1_air,ptfsn1_exp=read_comp(obj_r2)
else:
	ptfstar1_name,ptfstar1_air,ptfstar1_exp=read_comp(obj_r2)

#Blue parts	
out_b1=issn(obj_b1)+'2.fits'

if 'sn' in out_b1:
	ptfsn2_name,ptfsn2_air,ptfsn2_exp=read_comp(obj_b1)
else:
	ptfstar2_name,ptfstar2_air,ptfstar2_exp=read_comp(obj_b1)

out_b2=issn(obj_b2)+'2.fits'

if 'sn' in out_b2:
	ptfsn2_name,ptfsn2_air,ptfsn2_exp=read_comp(obj_b2)
else:
	ptfstar2_name,ptfstar2_air,ptfstar2_exp=read_comp(obj_b2)


print out_r1,out_r2,out_b1,out_b2

obj_r1='tgs'+obj_r1+'[SCI]'
obj_r2='tgs'+obj_r2+'[SCI]'
obj_b1='tgs'+obj_b1+'[SCI]'
obj_b2='tgs'+obj_b2+'[SCI]'

#Copy into new fits files (ptfsn1,ptfsn2,ptfstar1,ptfstar2)
iraf.imcopy(input=obj_r1,output=out_r1)
iraf.imcopy(input=obj_r2,output=out_r2)
iraf.imcopy(input=obj_b1,output=out_b1)
iraf.imcopy(input=obj_b2,output=out_b2)


#The new filenames
ptfsn1='ptfsn1.fits'
ptfsn2='ptfsn2.fits'
ptfstar1='ptfstar1.fits'
ptfstar2='ptfstar2.fits'

#Cosmic Cleaned Version
ptfsn1cr='ptfsn1.cr'
ptfsn2cr='ptfsn2.cr'

#Cosmic Rejection
iraf.lacos_spec(input=ptfsn1, output=ptfsn1cr, outmask='ptfsn1.crmask', gain=1.0,readn=3.85, xorder=9, yorder=3, sigclip=4.5, sigfrac=0.5, objlim=1.0, niter=3, verbose=yes, mode='al')

iraf.lacos_spec(input=ptfsn2, output=ptfsn2cr, outmask='ptfsn2.crmask', gain=1.0, readn=3.85, xorder=9, yorder=3, sigclip=4.5, sigfrac=0.5, objlim=1.0, niter=3, verbose=yes, mode='al')

#multispec
ptfsn1ms='ptfsn1.ms'
ptfsn2ms='ptfsn2.ms'
ptfstar1ms='ptfstar1.ms'
ptfstar2ms='ptfstar2.ms'

#Reduce with apall
iraf.apall(input=ptfsn1cr, out=ptfsn1ms, nfind=1, interactive=yes, find=yes, recenter=yes, resize=yes, edit=yes, trace=yes, fittrace=yes, t_nsum = 50, t_nlost = 50, extract=yes, extras=yes, review=yes, background='fit')

iraf.apall(input=ptfsn2cr, out=ptfsn2ms, nfind=1, interactive=yes, find=yes, recenter=yes, resize=yes, edit=yes, trace=yes, fittrace=yes, t_nsum = 50, t_nlost = 50, extract=yes, extras=yes, review=yes, background='fit')

iraf.apall(input=ptfstar1, out=ptfstar1ms)
iraf.apall(input=ptfstar2, out=ptfstar2ms)

#log?
l.write('extracted=yes')


print '\n'
print 'Done with reduction, now calibrate'
print '\n'

stdstar1='stdstar1'
stdstar2='stdstar2'

#Identify regions where we want to fit sensitivity function

iraf.standard(input=ptfstar1ms, output=stdstar1, star_name=ptfstar1_name, airmass=ptfstar1_air, exptime=ptfstar1_exp, answer=yes, caldir=default_caldir, extinction='mk_extinct.txt', observatory='Keck')

iraf.standard(input=ptfstar2ms, output=stdstar2, star_name=ptfstar2_name, airmass=ptfstar2_air, exptime=ptfstar2_exp, answer=yes, caldir=default_caldir, extinction='mk_extinct.txt', observatory='Keck')

sensstar1='sensstar1'
sensstar2='sensstar2'

#Fit Sensitivity Function
iraf.sensfunc(standards=stdstar1, sensitivity=sensstar1, observatory='Keck', order=4)
iraf.sensfunc(standards=stdstar2, sensitivity=sensstar2, observatory='Keck', order=4) 

ptfsn1f='ptfsn1.f'
ptfsn2f='ptfsn2.f'
sensstar1='sensstar1.fits'
sensstar2='sensstar2.fits'

#Calibrate Spectra
print '\n provide the following manually: \n'
print ptfsn1_name,ptfsn1_air,ptfsn1_exp
print '\n'

iraf.onedspec.calibrate(input=ptfsn1ms, out=ptfsn1f, extinction='mk_extinct.txt', observatory='Keck', ignoreaps=yes, sensitivity='sensstar1.fits')

print '\n provide the following manually: \n'
print ptfsn2_name,ptfsn2_air,ptfsn2_exp
print '\n'

iraf.onedspec.calibrate(input=ptfsn2ms, out=ptfsn2f, extinction='mk_extinct.txt', observatory='Keck', ignoreaps=yes, sensitivity='sensstar2.fits')

inputname=ptfsn1f+','+ptfsn2f
#Produce combined spectrum
iraf.scombine(input=inputname, out='finalspectrum.ms', reject='avsigclip', scale='median', sample='5500:6500')

l.write('calib=yes')

iraf.splot(images='finalspectrum.ms.fits')

iraf.wspectext(input='finalspectrum.ms.fits[*,1]',output=ptfsn1_name+'.ascii',header='NO')
iraf.wspectext(input='ptfsn1.f.fits[*,1]',output=ptfsn1_name+'_r400.ascii',header='NO')
iraf.wspectext(input='ptfsn2.f.fits[*,1]',output=ptfsn1_name+'_b600.ascii',header='NO')

l.write('finalname='+ptfsn1_name+'.ascii')

l.close()