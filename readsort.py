import astropy.io.fits as fits
import glob
import subprocess


#cmd='!gethead N*.fits -x 0 UT OBJECT GRATING EXPTIME | grep R400 > R400+_G5305.lst'
#subprocess.Popen(cmd).wait()
#cmd='!gethead N*.fits -x 0 UT OBJECT GRATING EXPTIME | grep B600 > B600+_G5307.lst'
#subprocess.Popen(cmd).wait()

files= glob.glob('N*.fits')

#for i in xrange(len(files)):
#	files[i]


f=open('R400+_G5305.lst','w')
s=open('B600+_G5307.lst','w')
w=1
i=1
j=1
k=1



for name in files:

    hdulist= fits.open(name)
    hdulist.close()
    if hdulist[0].header['grating']=='R400+_G5305':
        print>>f, name, hdulist[0].header['object']
        if 'flat' in hdulist[0].header['object']:
            x=open('flatr400_'+str(i)+'.txt','w')
            print>>x, name
            x.close()
            i+=1

        elif hdulist[0].header['OBSTYPE']=='OBJECT':
            x=open(str(hdulist[0].header['object'])+'_r400_'+str(j)+'.txt','w')
	    b=open('object_r400_'+str(w)+'.txt','w')
            print>>x,name
	    print>>b,name
	    b.close()
            j+=1
	    w+=1
            x.close()
        elif 'Ar' in hdulist[0].header['object']:
            x=open(str(hdulist[0].header['object'])+'_r400_'+str(k)+'.txt','w')
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
    if hdulist[0].header['grating']=='B600+_G5307':
        print>>s, name, hdulist[0].header['object']
        if 'flat' in hdulist[0].header['object']:
            x=open('flatb600_'+str(i)+'.txt','w')
            print>>x, name
            x.close()
            i+=1

        elif hdulist[0].header['OBSTYPE']=='OBJECT':
            x=open(str(hdulist[0].header['object'])+'_b600_'+str(j)+'.txt','w')
	    b=open('object_b600_'+str(w)+'.txt','w')
            print>>x,name
 	    print>>b,name
            j+=1
	    w+=1
            x.close()
	    b.close()
        elif 'Ar' in hdulist[0].header['object']:
            x=open(str(hdulist[0].header['object'])+'_b600_'+str(k)+'.txt','w')
            print>>x,name
            k+=1
            x.close()
f.close()
s.close()