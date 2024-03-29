import astropy.io.fits as fits
import glob
import subprocess

#!gethead N*.fits -x 0 UT OBJECT GRATING EXPTIME | grep B600 > B600+_G5307.lst
files= glob.glob('N*.fits')
f=open('B600+_G5307_520.lst','w')
s=open('B600+_G5307_525.lst','w')
w=1
i=1
j=1
k=1
for name in files:

    hdulist= fits.open(name)
    hdulist.close()
    if hdulist[0].header['CENTWAVE']==520:
        print>>f, name, hdulist[0].header['object']
        if 'flat' in hdulist[0].header['object']:
            x=open('flatb600_520_'+str(i)+'.txt','w')
            print>>x, name
            x.close()
            i+=1

        elif hdulist[0].header['OBSTYPE']=='OBJECT':
            x=open(str(hdulist[0].header['object'])+'_b600_520_'+str(j)+'.txt','w')
	    b=open('object_b600_520_'+str(w)+'.txt','w')
            print>>x,name
	    print>>b,name
	    b.close()
            j+=1
	    w+=1
            x.close()
        elif 'Ar' in hdulist[0].header['object']:
            x=open(str(hdulist[0].header['object'])+'_b600_520_'+str(k)+'.txt','w')
            print>>x,name
            k+=1
            x.close()
i=1
j=1
k=1
w=1
for name in files:

    hdulist= fits.open(name)
    hdulist.close()
    if hdulist[0].header['CENTWAVE']==525:
        print>>s, name, hdulist[0].header['object']
        if 'flat' in hdulist[0].header['object']:
            x=open('flatb600_525_'+str(i)+'.txt','w')
            print>>x, name
            x.close()
            i+=1

        elif hdulist[0].header['OBSTYPE']=='OBJECT':
            x=open(str(hdulist[0].header['object'])+'_b600_525_'+str(j)+'.txt','w')
	    b=open('object_b600_525_'+str(w)+'.txt','w')
            print>>x,name
 	    print>>b,name
            j+=1
	    w+=1
            x.close()
	    b.close()
        elif 'Ar' in hdulist[0].header['object']:
            x=open(str(hdulist[0].header['object'])+'_b600_525_'+str(k)+'.txt','w')
            print>>x,name
            k+=1
            x.close()
f.close()
s.close()
