% procedure to sort files for gemini reduction, reduce and calibrate spectra

set(0,'defaulttextinterpreter','latex')

cleanup='no' % clean up junk/temp files
lacosm='no' % enable/disable cosmic ray removal

if strcmp(cleanup,'yes')
	!rm tg*.fits
	!rm g*.fits
	!rm *.txt
	!rm normflat*
	!rm tmp*
	!rm OBJECT*.fits
	!rm r400*.fits
	!rm b600*.fits
	!rm lacos*
	!rm gsmask*
end



% collect relevant .fits header information from the cal. and science frames
!gethead N*.fits -x 0 UT OBJECT GRATING OBSTYPE > fitsheader.txt

   sn_out=fopen(strcat('fitsheader.txt'));
   temp_head=textscan(sn_out,'%s%s%s%s%s%s');
   fclose(sn_out);
   
   
   
% sort flats
   R400_fits=(strcmp(temp_head{5},'R400+_G5305')); 
   R_400_flat_t=(strcmp(temp_head{6},'FLAT')); 
   R_400_flat=temp_head{1}(boolean(R400_fits.*R_400_flat_t));
   
for n=1:length(R_400_flat)
	dlmwrite(strcat('flat_r400_',num2str(n),'.txt'),R_400_flat(n),'delimiter','')
end
   
   
   R600_fits=(strcmp(temp_head{5},'B600+_G5307'));  
   B_600_flat_t=(strcmp(temp_head{6},'FLAT'));  
   B_600_flat=temp_head{1}(boolean(R600_fits.*B_600_flat_t));
   
for n=1:length(B_600_flat)
	dlmwrite(strcat('flat_b600_',num2str(n),'.txt'),B_600_flat(n),'delimiter','')
end


% sort science frames (both std star and transient)
   R_400_object_t=(strncmp(temp_head{6},'OBJECT',4));
   R_400_obj=temp_head{1}(boolean(R400_fits.*R_400_object_t));   
   temp_h_sci400=temp_head{6}(boolean(R400_fits.*R_400_object_t));
   
for n=1:length(R_400_obj)
	dlmwrite(strcat(char(temp_h_sci400(n)),char(num2str(n)),'_400.txt'),R_400_obj(n),'delimiter','')
	dlmwrite(strcat(char(temp_h_sci400(n)),char(num2str(n)),'_400_ap.txt'),strcat(R_400_obj(n),'[SCI]'),'delimiter','')
end
   
   
   R_600_object_t=(strncmp(temp_head{6},'OBJECT',4));   
   R_600_obj=temp_head{1}(boolean(R600_fits.*R_600_object_t));
   temp_h_sci600=temp_head{6}(boolean(R600_fits.*R_600_object_t));
   
for n=1:length(R_600_obj)
	dlmwrite(strcat(char(temp_h_sci600(n)),char(num2str(n)),'_600.txt'),R_600_obj(n),'delimiter','')
	dlmwrite(strcat(char(temp_h_sci600(n)),char(num2str(n)),'_600_ap.txt'),strcat(R_600_obj(n),'[SCI]'),'delimiter','')
end
   
   
% arcs
   R_400_object_t=(strncmp(temp_head{6},'ARC',4));
   R_400_obj=temp_head{1}(boolean(R400_fits.*R_400_object_t));
   temp_h_arc400=temp_head{4}(boolean(R400_fits.*R_400_object_t));
   
for n=1:length(R_400_obj)
	dlmwrite(strcat(char(temp_h_arc400(n)),char(num2str(n)),'_400.txt'),R_400_obj(n),'delimiter','')
end
   
   
   R_600_object_t=(strncmp(temp_head{6},'ARC',4));
   R_600_obj=temp_head{1}(boolean(R600_fits.*R_600_object_t));   
   temp_h_arc600=temp_head{4}(boolean(R600_fits.*R_600_object_t));

for n=1:length(R_600_obj)
	dlmwrite(strcat(char(temp_h_arc600(n)),char(num2str(n)),'_600.txt'),R_600_obj(n),'delimiter','')
end
   
   
 
   
   %										 %
   % write and run pyraf gemini/gmos scripts %
   %										 %
   
   
%make norm. flats
for n=1:length(R_400_flat)
	   
	p1=sprintf('from pyraf import iraf\nfrom pyraf.iraf import gemini\nfrom pyraf.iraf import gmos\n\n');
	p2=sprintf(strcat('iraf.gemini.gsflat(inflats="@',strcat('flat_r400_',num2str(n),'.txt'),'",specflat="',strcat('normflat_r400_',num2str(n),'.fits'),'",fl_bias="no",fl_inter="yes",fl_detec="yes",fl_seprows="no")'));
	dlmwrite('gemini_reduce.py',{p1 p2},'delimiter','');
		dlmwrite('gemini_reduce1.py',{p1 p2},'delimiter','');

   
	cmd=sprintf('/bin/bash --login -c ''ur_setup; ipython gemini_reduce.py''');
    %[r,s]=
	[r,s]=system(cmd);
end
   
for n=1:length(B_600_flat)
	p1=sprintf('from pyraf import iraf\nfrom pyraf.iraf import gemini\nfrom pyraf.iraf import gmos\n\n');
	p2=sprintf(strcat('iraf.gemini.gsflat(inflats="@',strcat('flat_b600_',num2str(n),'.txt'),'",specflat="',strcat('normflat_b600_',num2str(n),'.fits'),'",fl_bias="no",fl_inter="yes",fl_detec="yes",fl_seprows="no")'));
	dlmwrite('gemini_reduce.py',{p1 p2},'delimiter','');
   		dlmwrite('gemini_reduce2.py',{p1 p2},'delimiter','');

    cmd=sprintf('/bin/bash --login -c ''ur_setup; ipython gemini_reduce.py''');
	[r,s]=system(cmd);
end
   
% sort object frames to correct flat fields
% this needs to be fixed to match flat field and object numbers by matching by closest time of obs (UT fitsheader)   
for n=1:length(R_400_flat)
	normflat400{n}=strcat('normflat_r400_',num2str(n),'.fits');
	normflat_arc400{n}=strcat('normflat_r400_',num2str(n),'.fits');
end
   
for n=1:length(B_600_flat)
	normflat600{n}=strcat('normflat_b600_',num2str(n),'.fits');
	normflat_arc600{n}=strcat('normflat_b600_',num2str(n),'.fits');
end
   
for n=1:length(temp_h_sci400)
	if n > size(normflat400,2)
		normflat400{n}=normflat400{n-1};
	end
end
   
for n=1:length(temp_h_arc400)
	if n > size(normflat400,2)
		normflat_arc400{n}=normflat_arc400{n-1};
	end
end
   
for n=1:length(temp_h_sci600)
	if n > size(normflat600,2)
		normflat600{n}=normflat600{n-1};
	end
end
   
for n=1:length(temp_h_arc600)
	if n > size(normflat600,2)
		normflat_arc600{n}=normflat_arc600{n-1};
	end
end
   
   
%apply norm. flat to object frames (arc and science frames)
for n=1:length(temp_h_sci400)
	p1=sprintf('from pyraf import iraf\nfrom pyraf.iraf import gemini\nfrom pyraf.iraf import gmos\n\n');
	p2=sprintf(strcat('iraf.gemini.gsreduce(inimages="@',strcat(char(temp_h_sci400(n)),char(num2str(n)),'_400.txt'),'",flatim="',char(normflat400{n}),'",fl_bias="no")'));
	dlmwrite('gemini_reduce.py',{p1 p2},'delimiter','');
			dlmwrite('gemini_reduce3.py',{p1 p2},'delimiter','');

   
	cmd=sprintf('/bin/bash --login -c ''ur_setup; ipython gemini_reduce.py''');
    system(cmd)
end
   
for n=1:length(temp_h_sci600)
	p1=sprintf('from pyraf import iraf\nfrom pyraf.iraf import gemini\nfrom pyraf.iraf import gmos\n\n');
	p2=sprintf(strcat('iraf.gemini.gsreduce(inimages="@',strcat(char(temp_h_sci600(n)),char(num2str(n)),'_600.txt'),'",flatim="',char(normflat600{n}),'",fl_bias="no")')); 
	dlmwrite('gemini_reduce.py',{p1 p2},'delimiter','');
			dlmwrite('gemini_reduce4.py',{p1 p2},'delimiter','');

   
	cmd=sprintf('/bin/bash --login -c ''ur_setup; ipython gemini_reduce.py''');
    system(cmd)
end
	
	
for n=1:length(temp_h_arc400)
	p1=sprintf('from pyraf import iraf\nfrom pyraf.iraf import gemini\nfrom pyraf.iraf import gmos\n\n');
	p2=sprintf(strcat('iraf.gemini.gsreduce(inimages="@',strcat(char(temp_h_arc400(n)),char(num2str(n)),'_400.txt'),'",flatim="',char(normflat_arc400{n}),'",fl_bias="no")')); 
	dlmwrite('gemini_reduce.py',{p1 p2},'delimiter','');
				dlmwrite('gemini_reduce5.py',{p1 p2},'delimiter','');

   
	cmd=sprintf('/bin/bash --login -c ''ur_setup; ipython gemini_reduce.py''');
    system(cmd)
end
   
for n=1:length(temp_h_arc600)
	p1=sprintf('from pyraf import iraf\nfrom pyraf.iraf import gemini\nfrom pyraf.iraf import gmos\n\n');
	p2=sprintf(strcat('iraf.gemini.gsreduce(inimages="@',strcat(char(temp_h_arc600(n)),char(num2str(n)),'_600.txt'),'",flatim="',char(normflat_arc600{n}),'",fl_bias="no")'));
	dlmwrite('gemini_reduce.py',{p1 p2},'delimiter','');
				dlmwrite('gemini_reduce6.py',{p1 p2},'delimiter','');

   
	cmd=sprintf('/bin/bash --login -c ''ur_setup; ipython gemini_reduce.py''');
    system(cmd)
end
   
   
% solve for the wavelength solution 
for n=1:length(temp_h_arc400)
	p1=sprintf('from pyraf import iraf\nfrom pyraf.iraf import gemini\nfrom pyraf.iraf import gmos\n\n');
	p2=sprintf(strcat('iraf.gemini.gswavelength(inimages="gs//@',strcat(char(temp_h_arc400(n)),char(num2str(n)),'_400.txt'),'")'));   
	dlmwrite('gemini_reduce.py',{p1 p2},'delimiter','');
				dlmwrite('gemini_reduce7.py',{p1 p2},'delimiter','');

   
	cmd=sprintf('/bin/bash --login -c ''ur_setup; ipython gemini_reduce.py''');
	%f = warndlg('Push enter several times in the matlab terminal, if nothing happens.', 'Warn Dialog');
    system(cmd)
	%close(f)
end
   
for n=1:length(temp_h_arc600)
	p1=sprintf('from pyraf import iraf\nfrom pyraf.iraf import gemini\nfrom pyraf.iraf import gmos\n\n');
	p2=sprintf(strcat('iraf.gemini.gswavelength(inimages="gs//@',strcat(char(temp_h_arc600(n)),char(num2str(n)),'_600.txt'),'")'));
	dlmwrite('gemini_reduce.py',{p1 p2},'delimiter','');
					dlmwrite('gemini_reduce8.py',{p1 p2},'delimiter','');

   
	%f = warndlg('Push enter several times in the matlab terminal, if nothing happens.', 'Warn Dialog');
	cmd=sprintf('/bin/bash --login -c ''ur_setup; ipython gemini_reduce.py''');
    system(cmd)
	%close(f)
end

% match arcs to sci frames
% this needs to be fixed to match flat field and object (numbers) by matching closest times (UT fitsheader)
for n=1:length(temp_h_arc400)
	arc_match_400{n}=strcat(char(temp_h_arc400(n)),char(num2str(n)),'_400.txt');
end

for n=1:length(temp_h_sci400)
	if n > size(arc_match_400,2)
		arc_match_400{n}=arc_match_400{n-1};
	end
end

for n=1:length(temp_h_arc600)
	arc_match_600{n}=strcat(char(temp_h_arc600(n)),char(num2str(n)),'_600.txt');
end

for n=1:length(temp_h_sci600)
	if n > size(arc_match_600,2)
		arc_match_600{n}=arc_match_600{n-1};
	end
end
   
% apply the wavelength solutions
for n=1:length(temp_h_sci400)
	p1=sprintf('from pyraf import iraf\nfrom pyraf.iraf import gemini\nfrom pyraf.iraf import gmos\n\n');
	p2=sprintf(strcat('iraf.gemini.gstransform(inimages="gs//@',strcat(char(temp_h_sci400(n)),char(num2str(n)),'_400.txt'),'",wavtran="gs//@',arc_match_400{n},'")'));   
	dlmwrite('gemini_reduce.py',{p1 p2},'delimiter','');
   					dlmwrite('gemini_reduce9.py',{p1 p2},'delimiter','');

	cmd=sprintf('/bin/bash --login -c ''ur_setup; ipython gemini_reduce.py''');
    system(cmd)
end

for n=1:length(temp_h_sci600)
	p1=sprintf('from pyraf import iraf\nfrom pyraf.iraf import gemini\nfrom pyraf.iraf import gmos\n\n');
	p2=sprintf(strcat('iraf.gemini.gstransform(inimages="gs//@',strcat(char(temp_h_sci600(n)),char(num2str(n)),'_600.txt'),'",wavtran="gs//@',arc_match_600{n},'")'));   
	dlmwrite('gemini_reduce.py',{p1 p2},'delimiter','');
						dlmwrite('gemini_reduce10.py',{p1 p2},'delimiter','');

   
	cmd=sprintf('/bin/bash --login -c ''ur_setup; ipython gemini_reduce.py''');
    system(cmd)
end



% collect images in which to extract traces
   R_400_object_t=(strncmp(temp_head{6},'OBJECT',4));   
   R_400_obj=temp_head{1}(boolean(R400_fits.*R_400_object_t));  
   temp_h_sci400=temp_head{6}(boolean(R400_fits.*R_400_object_t));
   temp_h_sciout400=strcat('r400_',temp_head{4}(boolean(R400_fits.*R_400_object_t)));

   R_600_object_t=(strncmp(temp_head{6},'OBJECT',4));  
   R_600_obj=temp_head{1}(boolean(R600_fits.*R_600_object_t));  
   temp_h_sci600=temp_head{6}(boolean(R600_fits.*R_600_object_t));
   temp_h_sciout600=strcat('b600_',temp_head{4}(boolean(R600_fits.*R_600_object_t)));
   
if strcmp(lacosm,'yes')  
% cosmic ray removal enabled

for n=1:length(temp_h_sci400)
	p1=sprintf('from pyraf import iraf\nfrom pyraf.iraf import gemini\nfrom pyraf.iraf import gmos\n\n');
	p2=sprintf(strcat('iraf.lacos_spec(input="gs//@',strcat(char(temp_h_sci400(n)),char(num2str(n)),'_400_ap.txt'),'",output="lacos_',strcat(char(temp_h_sci400(n)),char(num2str(n)),'_400.fits'),'",outmask="gsmask//@',strcat(char(temp_h_sci400(n)),char(num2str(n)),'_400.txt'),'",gain="1.0",readn="3.85",xorder="9",yorder="3",sigclip="4.5",sigfrac="0.5",objlim="1.0",niter="3")'));
	dlmwrite('gemini_reduce.py',{p1 p2},'delimiter','')
						dlmwrite('gemini_reduce11.py',{p1 p2},'delimiter','');

	
   	%f = warndlg('Push enter several times in the matlab terminal, if nothing happens.', 'Warn Dialog');
	
	cmd=sprintf('/bin/bash --login -c ''ur_setup; ipython gemini_reduce.py''');
    system(cmd)
end

for n=1:length(temp_h_sci600)
	p1=sprintf('from pyraf import iraf\nfrom pyraf.iraf import gemini\nfrom pyraf.iraf import gmos\n\n');
	p2=sprintf(strcat('iraf.lacos_spec(input="gs//@',strcat(char(temp_h_sci600(n)),char(num2str(n)),'_600_ap.txt'),'",output="lacos_',strcat(char(temp_h_sci600(n)),char(num2str(n)),'_600.fits'),'",outmask="gsmask//@',strcat(char(temp_h_sci600(n)),char(num2str(n)),'_600.txt'),'",gain="1.0",readn="3.85",xorder="9",yorder="3",sigclip="4.5",sigfrac="0.5",objlim="1.0",niter="3")'));
	dlmwrite('gemini_reduce.py',{p1 p2},'delimiter','')
						dlmwrite('gemini_reduce12.py',{p1 p2},'delimiter','');

	
   	%f = warndlg('Push enter several times in the matlab terminal, if nothing happens.', 'Warn Dialog');
	
	cmd=sprintf('/bin/bash --login -c ''ur_setup; ipython gemini_reduce.py''');
    system(cmd)
end

% extract traces in image frames  
for n=1:length(temp_h_sci400)
	p1=sprintf('from pyraf import iraf\nfrom pyraf.iraf import gemini\nfrom pyraf.iraf import gmos\n\n');
	p2=sprintf(strcat('iraf.apall(input="lacos_',strcat(char(temp_h_sci400(n)),char(num2str(n)),'_400.fits'),'",out="',char(temp_h_sciout400{n}),num2str(n),'.fits",interactive="yes",find="yes",recenter="yes",resize="yes",edit="yes",trace="yes",fittrace="yes",extract="yes",extras="yes",review="yes",background="fit")'));
	dlmwrite('gemini_reduce.py',{p1 p2},'delimiter','')
						dlmwrite('gemini_reduce13.py',{p1 p2},'delimiter','');

	
   	%f = warndlg('Push enter several times in the matlab terminal, if nothing happens.', 'Warn Dialog');
	
	cmd=sprintf('/bin/bash --login -c ''ur_setup; ipython gemini_reduce.py''');
    system(cmd)
    %close(f)
	
	% write 1D extracted science and sky spectra
	p3=sprintf(strcat('iraf.wspectext(input="',strcat(char(temp_h_sciout400{n}),num2str(n),'.fits[*,1,1]'),'",output="',strcat(char(temp_h_sciout400{n}),num2str(n),'.txt'),'")\n'));
	p4=sprintf(strcat('iraf.wspectext(input="',strcat(char(temp_h_sciout400{n}),num2str(n),'.fits[*,1,2]'),'",output="',strcat(char(temp_h_sciout400{n}),num2str(n),'_sky.txt'),'")\n'));
	dlmwrite('gemini_reduce.py',{p1 p3 p4},'delimiter','');
		dlmwrite('gemini_reduce14.py',{p1 p3 p4},'delimiter','');

	
  	cmd=sprintf('/bin/bash --login -c ''ur_setup; ipython gemini_reduce.py''');
    system(cmd)
end

  
for n=1:length(temp_h_sci600)
	p1=sprintf('from pyraf import iraf\nfrom pyraf.iraf import gemini\nfrom pyraf.iraf import gmos\n\n');
	p2=sprintf(strcat('iraf.apall(input="gs//@',strcat(char(temp_h_sci600(n)),char(num2str(n)),'_600_ap.txt'),'",out="',char(temp_h_sciout600{n}),num2str(n),'.fits",nfind="1",interactive="yes",find="yes",recenter="yes",resize="yes",edit="yes",trace="yes",fittrace="yes",extract="yes",extras="yes",review="yes",background="fit")'));
	dlmwrite('gemini_reduce.py',{p1 p2},'delimiter','')
							dlmwrite('gemini_reduce15.py',{p1 p2},'delimiter','');

	
   	%f = warndlg('Push enter several times in the matlab terminal, if nothing happens.', 'Warn Dialog');
	
	cmd=sprintf('/bin/bash --login -c ''ur_setup; ipython gemini_reduce.py''');
    system(cmd)
	%close(f)

	% write 1D extracted science and sky spectra
	p3=sprintf(strcat('iraf.wspectext(input="',strcat(char(temp_h_sciout600{n}),num2str(n),'.fits[*,1,1]'),'",output="',strcat(char(temp_h_sciout600{n}),num2str(n),'.txt'),'")\n'));
    p4=sprintf(strcat('iraf.wspectext(input="',strcat(char(temp_h_sciout600{n}),num2str(n),'.fits[*,1,2]'),'",output="',strcat(char(temp_h_sciout600{n}),num2str(n),'_sky.txt'),'")\n'));
	dlmwrite('gemini_reduce.py',{p1 p3 p4},'delimiter','')
		dlmwrite('gemini_reduce16.py',{p1 p3 p4},'delimiter','')

	
	cmd=sprintf('/bin/bash --login -c ''ur_setup; ipython gemini_reduce.py''');
    system(cmd)
end

end


if strcmp(lacosm,'no') 
% cosmic ray removal disabled

% extract traces in image frames  
for n=1:length(temp_h_sci400)
	p1=sprintf('from pyraf import iraf\nfrom pyraf.iraf import gemini\nfrom pyraf.iraf import gmos\n\n');
	p2=sprintf(strcat('iraf.apall(input="gs//@',strcat(char(temp_h_sci400(n)),char(num2str(n)),'_400_ap.txt'),'",out="',char(temp_h_sciout400{n}),num2str(n),'.fits",interactive="yes",find="yes",recenter="yes",resize="yes",edit="yes",trace="yes",fittrace="yes",extract="yes",extras="yes",review="yes",background="fit")'));
	dlmwrite('gemini_reduce.py',{p1 p2},'delimiter','')
								dlmwrite('gemini_reduce17.py',{p1 p2},'delimiter','');

	
   	%f = warndlg('Push enter several times in the matlab terminal, if nothing happens.', 'Warn Dialog');
	
	cmd=sprintf('/bin/bash --login -c ''ur_setup; ipython gemini_reduce.py''');
    system(cmd)
    %close(f)
	
	% write 1D extracted science and sky spectra
	p3=sprintf(strcat('iraf.wspectext(input="',strcat(char(temp_h_sciout400{n}),num2str(n),'.fits[*,1,1]'),'",output="',strcat(char(temp_h_sciout400{n}),num2str(n),'.txt'),'")\n'));
	p4=sprintf(strcat('iraf.wspectext(input="',strcat(char(temp_h_sciout400{n}),num2str(n),'.fits[*,1,2]'),'",output="',strcat(char(temp_h_sciout400{n}),num2str(n),'_sky.txt'),'")\n'));
	dlmwrite('gemini_reduce.py',{p1 p3 p4},'delimiter','');
		dlmwrite('gemini_reduce18.py',{p1 p3 p4},'delimiter','');

	
  	cmd=sprintf('/bin/bash --login -c ''ur_setup; ipython gemini_reduce.py''');
    system(cmd)
end

  
for n=1:length(temp_h_sci600)
	p1=sprintf('from pyraf import iraf\nfrom pyraf.iraf import gemini\nfrom pyraf.iraf import gmos\n\n');
	p2=sprintf(strcat('iraf.apall(input="gs//@',strcat(char(temp_h_sci600(n)),char(num2str(n)),'_600_ap.txt'),'",out="',char(temp_h_sciout600{n}),num2str(n),'.fits",nfind="1",interactive="yes",find="yes",recenter="yes",resize="yes",edit="yes",trace="yes",fittrace="yes",extract="yes",extras="yes",review="yes",background="fit")'));
	dlmwrite('gemini_reduce.py',{p1 p2},'delimiter','')
									dlmwrite('gemini_reduce19.py',{p1 p2},'delimiter','');

	
   	%f = warndlg('Push enter several times in the matlab terminal, if nothing happens.', 'Warn Dialog');
	
	cmd=sprintf('/bin/bash --login -c ''ur_setup; ipython gemini_reduce.py''');
    system(cmd)
	%close(f)

	% write 1D extracted science and sky spectra
	p3=sprintf(strcat('iraf.wspectext(input="',strcat(char(temp_h_sciout600{n}),num2str(n),'.fits[*,1,1]'),'",output="',strcat(char(temp_h_sciout600{n}),num2str(n),'.txt'),'")\n'));
    p4=sprintf(strcat('iraf.wspectext(input="',strcat(char(temp_h_sciout600{n}),num2str(n),'.fits[*,1,2]'),'",output="',strcat(char(temp_h_sciout600{n}),num2str(n),'_sky.txt'),'")\n'));
	dlmwrite('gemini_reduce.py',{p1 p3 p4},'delimiter','')
		dlmwrite('gemini_reduce20.py',{p1 p3 p4},'delimiter','')

	
	cmd=sprintf('/bin/bash --login -c ''ur_setup; ipython gemini_reduce.py''');
    system(cmd)
end

end






%										   %
% flux calibration and telluric correction %
%									       %


% find correct standard star spectrum for the observation
for n=1:size(temp_h_sciout400,1)
	if strncmp(temp_h_sciout400{n}(6:end),'iPTF',4)==0
		%standard_n=n;
		standard=search_specphot_stand(char(temp_h_sciout400{n}(6:end)),'substr');
		plotdata_std=dlmread(strcat(char(temp_h_sciout400{n}),num2str(n),'.txt'));	
		break % <-- now we just pick the first std. star observation of the run, should show all of them and allow for manual choosing if two different standards were observed
	end
end


wl400=[5500:0.2:9200];
wl600=[3780:0.2:6500];

% plot raw extracted spectra of sci. target and std

standard_n400=[];
for n=1:size(temp_h_sciout400,1)
	if strncmp(temp_h_sciout400{n}(6:end),'iPTF',4)==0
		standard_n400=[standard_n400 n];
	end
end

standard_n600=[];
for n=1:size(temp_h_sciout600,1)
	if strncmp(temp_h_sciout600{n}(6:end),'iPTF',4)==0
		standard_n600=[standard_n600 n];
	end
end


k=1;
n=1:size(temp_h_sciout600,1);
match=setxor(n,standard_n600);

for n=match
			plotdata=dlmread(strcat(char(temp_h_sciout600{n}),num2str(n),'.txt'));
	
			figure(1)
			subplot((size(temp_h_sciout600,1)+size(temp_h_sciout400,1))-(length(standard_n600)+length(standard_n400)),1,k)
			plot(plotdata(:,1),plotdata(:,2),'k-','linewidth',1)
			xlabel('Wavelength [Angstr]')
			ylabel('Flux (arb. units)')
			box on
			k=k+1;
end

k=length(standard_n600)+1;
n=1:size(temp_h_sciout400,1);
match=setxor(n,standard_n400);

for n=match
			plotdata=dlmread(strcat(char(temp_h_sciout400{n}),num2str(n),'.txt'));
			subplot((size(temp_h_sciout600,1)+size(temp_h_sciout400,1))-(length(standard_n600)+length(standard_n400)),1,k)
			plot(plotdata(:,1),plotdata(:,2),'k-','linewidth',1)
			xlabel('Wavelength [Angstr]')
			ylabel('Flux (arb. units)')
			box on
			k=k+1;
end

% prompt to choose best std. star observation, in case several exposures were made using different exp. times
figure(2)
		
	subplot(max([length(standard_n600) length(standard_n400)])+1,2,1)
	plot(standard.Spec(:,1),standard.Spec(:,2),'r-','linewidth',1)
	xlabel('Wavelength [Angstr]')
	ylabel('Flux (erg/cm2/s)')
	title(strcat('Grism 400, Spectrophot, standard',{','},{' '},char(standard.Name(1))),'fontsize',18)	
	xlim([4900 9400])
		
	subplot(max([length(standard_n600) length(standard_n400)])+1,2,2)
	plot(standard.Spec(:,1),standard.Spec(:,2),'r-','linewidth',1)
	xlabel('Wavelength [Angstr]')
	ylabel('Flux (erg/cm2/s)')
	title(strcat('Grism 600, Spectrophot, standard',{','},{' '},char(standard.Name(1))),'fontsize',18)	
	xlim([3700 6700])
		
		
k=1;
n=1:size(temp_h_sciout400,1);
match=intersect(n,standard_n400);

for n=match
	
	plotdata=dlmread(strcat(char(temp_h_sciout400{n}),num2str(n),'.txt'));
	temp400=textread(strcat(char(temp_h_sci400{n}),num2str(n),'_400.txt'),'%s');
	AirMass=cell2mat(get_fits_keyword(temp400{1},{'AIRMASS'},[]));
	
	EXPTIME=cell2mat(get_fits_keyword(temp400{1},{'EXPTIME'},[]));
	plotdata(:,2)=plotdata(:,2)./EXPTIME;
	
	plotdata=atmospheric_ext(plotdata,AirMass,'mkoextinct.dat','linear'); % atm. ext correction
	
	subplot(max([length(standard_n600) length(standard_n400)])+1,2,2.*k+1)		
	wlt{n}=plotdata(:,1);
	std_spec{n}=plotdata(:,2);
	plot(wlt{n},std_spec{n},'k-','linewidth',1)
	xlabel('Wavelength [Angstr]')
	ylabel('Flux (arb. units)')
	title(strcat('Standard star obs.n',{', '},'object nr',{' '},num2str(standard_n400(k)),{','},{' '},char(standard.Name(1))),'fontsize',18)		
	box on	
	xlim([4900 9400])
			
	k=k+1;
end

k=1;
n=1:size(temp_h_sciout600,1);
match=intersect(n,standard_n600);

for n=match
	plotdata=dlmread(strcat(char(temp_h_sciout600{n}),num2str(n),'.txt'));
	temp600=textread(strcat(char(temp_h_sci600{n}),num2str(n),'_600.txt'),'%s');
	AirMass=cell2mat(get_fits_keyword(temp600{1},{'AIRMASS'},[]));
	
	EXPTIME=cell2mat(get_fits_keyword(temp600{1},{'EXPTIME'},[]));
	plotdata(:,2)=plotdata(:,2)./EXPTIME;
			
	plotdata=atmospheric_ext(plotdata,AirMass,'mkoextinct.dat','linear'); % atm. ext correction
	
	subplot(max([length(standard_n600) length(standard_n400)])+1,2,2.*k+2)		
	wlt{n}=plotdata(:,1);
	std_spec{n}=plotdata(:,2);
	plot(wlt{n},std_spec{n},'k-','linewidth',1)
	xlabel('Wavelength [Angstr]')
	ylabel('Flux (arb. units)')
	title(strcat('Standard star obs.n',{', '},'object nr',{' '},num2str(standard_n600(k)),{','},{' '},char(standard.Name(1))),'fontsize',18)		
	box on	
	xlim([3700 6700])
			
	k=k+1;
end




if length(standard_n400)>1
x = inputdlg('Which standard star spectrum do you want to use for the r400 grating? input the oject number',...
	'Sample', [1 50]);
	standard_n400 = str2num(x{:});	
end


if length(standard_n600)>1
	x = inputdlg('Which standard star spectrum do you want to use for the r400 grating? input the oject number',...
	'Sample', [1 50]);
	standard_n600 = str2num(x{:});	
end


% find sensitivity function
for n=standard_n400;
	obsdata=dlmread(strcat(char(temp_h_sciout400{n}),num2str(n),'.txt'));
	temp400=textread(strcat(char(temp_h_sci400{n}),num2str(n),'_400.txt'),'%s');
	AirMass=cell2mat(get_fits_keyword(temp400{1},{'AIRMASS'},[]));
	
	EXPTIME=cell2mat(get_fits_keyword(temp400{1},{'EXPTIME'},[]));
	obsdata(:,2)=obsdata(:,2)./EXPTIME;
	
	obsdata=atmospheric_ext(obsdata,AirMass,'mkoextinct.dat','linear'); % atm. ext correction
	
	interpstandard=interp1(standard.Spec(:,1),standard.Spec(:,2),wl400,'poly1');
	interpobject=interp1(obsdata(:,1),obsdata(:,2),wl400,'poly1');
	
% % % 	[I interpobject_t]=clip_resid(interpobject,'Method','stdp','clip',[3 3]);
% % % 	interpstandard_t=interpstandard(I);
% % % 	wl400_t=wl400(I);

% % % 	calfactors=smooth(interpstandard,3)./smooth(interpobject,3);
% % % 	calfit=fit(wl400',calfactors,'poly3','normalize','on');
% % % 	calfactors=calfit(wl400);
	
	%exclude atm. abs lines, etc from sens. func. fit
	stdfit=fit(wl400',smooth(interpstandard,3),'poly7','normalize','on');
	objfit=fit(wl400',smooth(interpobject,3),'poly7','normalize','on');
	
	%exclude
		interp_s_hi=interpstandard'<1.2.*stdfit(wl400);
		interp_s_lo=interpstandard'>0.85.*stdfit(wl400);
		interp_s_i=boolean(interp_s_hi.*interp_s_lo);
	
		interp_o_hi=interpobject'<1.2.*objfit(wl400);
		interp_o_lo=interpobject'>0.85.*objfit(wl400);
		interp_o_i=boolean(interp_o_hi.*interp_o_lo);
	
		%interp_i=boolean(interp_o_i.*interp_s_i);
		
		% grow rejection by -10+20 units
		interp_o_i_temp=interp_o_i;
		for i=1:length(interp_o_i)
			if interp_o_i(i)==0 & i > 50 & i < length(interp_o_i)-150
				interp_o_i_temp(i-50:1:i+150)=0;
			end
		end
		
		% grow rejection by -10+20 units
		interp_s_i_temp=interp_s_i;
		for i=1:length(interp_s_i)
			if interp_s_i(i)==0 & i > 50 & i < length(interp_s_i)-150
				interp_s_i_temp(i-50:1:i+150)=0;
			end
		end
	
		interp_o_i=interp_o_i_temp;
		interp_s_i=interp_s_i_temp;
		
		stdfit=fit(wl400(interp_s_i)',smooth(interpstandard(interp_s_i),1),'poly7','normalize','on');
		
		% divide spectrum in 3 parts and fit polynomials

		wl_f=wl400(interp_o_i)';
		interpobject_f=interpobject(interp_o_i);
	
		

for m=1:999	
		figure(666)
		clf
		hold on
		plot(wl400(interp_o_i),interpobject(interp_o_i),'ko','markersize',2,'markerfacecolor',rgb('SkyBlue'),'color',rgb('SkyBlue'))
		plot(wl400(~interp_o_i),interpobject(~interp_o_i),'ko','markersize',2,'markerfacecolor',rgb('Salmon'),'color',rgb('Salmon'))
		plot(wl400,objfit(wl400),'--','color',rgb('Silver'),'linewidth',1.5)
		title('Sensitiviity function fits, and excluded regions, Observed Spectr. standard')
		box on
		
		x = inputdlg('Enter order of 3 polynomials to fit, e.g. 3 3 3:',...
		'Sample', [1 50]);
		orders = str2num(x{:});	
		
		'Select cutoff points for 3 polynomials to fit'
		[input_m,y] = ginput(2);
		input1=input_m(1);
		input2=input_m(2);
		

		
		wl_f1=wl_f(1:2000+round(mean(wl_f(wl_f<=input1.*1.01 & wl_f >=input1.*0.99))));
		wl_f2=wl_f(-2000+round(mean(wl_f(wl_f<=input1.*1.01 & wl_f >=input1.*0.99))):2000+round(mean(wl_f(wl_f<=input2.*1.01 & wl_f >=input2.*0.99))));
		wl_f3=wl_f(-2000+round(mean(wl_f(wl_f<=input2.*1.01 & wl_f >=input2.*0.99))):end);

		o_f1=smooth(interpobject_f(1:2000+round(mean(wl_f(wl_f<=input1.*1.01 & wl_f >=input1.*0.99)))),1);		
		o_f2=smooth(interpobject_f(-2000+round(mean(wl_f(wl_f<=input1.*1.01 & wl_f >=input1.*0.99))):2000+round(mean(wl_f(wl_f<=input2.*1.01 & wl_f >=input2.*0.99)))),1);
		o_f3=smooth(interpobject_f(-2000+round(mean(wl_f(wl_f<=input2.*1.01 & wl_f >=input2.*0.99))):end),1);
		
		objfit1=fit(wl_f1,o_f1,strcat('poly',num2str(orders(1))),'normalize','on');	
		objfit2=fit(wl_f2,o_f2,strcat('poly',num2str(orders(2))),'normalize','on');	
		objfit3=fit(wl_f3,o_f3,strcat('poly',num2str(orders(3))),'normalize','on');	

		objfit_p1=objfit1(wl400(1):0.2:round(input1));
		objfit_p2=objfit2(round(input1)+0.2:0.2:round(input2));
		objfit_p3=objfit3(round(input2)+0.2:0.2:wl400(end));
		
		stdfit_p=smooth(stdfit(wl400),1000);
		objfit_p=smooth([objfit_p1;objfit_p2;objfit_p3],1000);

		calfactors=stdfit_p./objfit_p;
		
		
		wl_tell=wl400(~interp_o_i)';
		interpobject_tell=interpobject(~interp_o_i);
		interpstandard_tell=interpstandard(~interp_o_i);
		calfactors_tell_temp=interpstandard_tell./interpobject_tell;
		calfactors_tell=calfactors.*1;
		calfactors_tell(~interp_o_i)=calfactors_tell_temp;
		
		% telluric correction
		calfactors=calfactors_tell;
	
		figure(666)
		hold on
		plot(wl400,objfit_p,'-','color',rgb('DarkOrchid'),'linewidth',1.5)
		title('Sensitiviity function fits, and excluded regions, Observed Spectr. standard')
		box on
	
		m=m+1;	
		%m = input('If fit is satisfactory type 0 and press enter. Otherwise press enter to redo fit.');
		x = inputdlg('If fit is satisfactory type 0 and press ok. Otherwise press ok to redo fit.',...
		'Sample', [1 50]);
		m = str2num(x{:});	
		if m==0
			break
		end
end

	figure(66)
	subplot(2,3,1)
	hold on
	plot(wl400(interp_s_i),interpstandard(interp_s_i),'ko','markersize',2,'markerfacecolor',rgb('SkyBlue'),'color',rgb('SkyBlue'))
	plot(wl400(~interp_s_i),interpstandard(~interp_s_i),'ko','markersize',2,'markerfacecolor',rgb('Salmon'),'color',rgb('Salmon'))
	plot(wl400,stdfit_p,'-','color',rgb('DarkOrchid'),'linewidth',1.5)
	plot(wl400,stdfit(wl400),'--','color',rgb('Silver'),'linewidth',1.5)

	title('Sensitiviity function fits, and excluded regions, Database Spectr. standard')
	box on

	subplot(2,3,2)
	hold on
	plot(wl400(interp_o_i),interpobject(interp_o_i),'ko','markersize',2,'markerfacecolor',rgb('SkyBlue'),'color',rgb('SkyBlue'))
	plot(wl400(~interp_o_i),interpobject(~interp_o_i),'ko','markersize',2,'markerfacecolor',rgb('Salmon'),'color',rgb('Salmon'))
	plot(wl400,objfit_p,'-','color',rgb('DarkOrchid'),'linewidth',1.5)
	plot(wl400,objfit(wl400),'--','color',rgb('Silver'),'linewidth',1.5)
	title('Sensitiviity function fits, and excluded regions, Observed Spectr. standard')
	box on
	
	subplot(2,3,3)
	plot(wl400,calfactors)
	title('Sensitivity function')
	box on

	% apply to science frame(s)
	for n=1:size(temp_h_sciout400,1)
		if strncmp(temp_h_sciout400{n}(6:end),'iPTF',4)==1
			obsdata=dlmread(strcat(char(temp_h_sciout400{n}),num2str(n),'.txt'));
			temp400=textread(strcat(char(temp_h_sci400{n}),num2str(n),'_400.txt'),'%s');
			AirMass=cell2mat(get_fits_keyword(temp400{1},{'AIRMASS'},[]));
			
			EXPTIME=cell2mat(get_fits_keyword(temp400{1},{'EXPTIME'},[]));
			obsdata(:,2)=obsdata(:,2)./EXPTIME;
			
			obsdata=atmospheric_ext(obsdata,AirMass,'mkoextinct.dat','linear'); % atm. ext correction
	
			interpobject=interp1(obsdata(:,1),obsdata(:,2),wl400,'poly1');

	figure(100+n)
	flux400=interpobject.*calfactors';
			plot(wl400,flux400,'k-','linewidth',1)
			xlabel('Wavelength [Angstr]')
			ylabel('Flux (erg/cm2/s)')
			title('Flux calibrated spectrum')
			box on
		end
	end
	
		% apply to std frame(s)
	for n=1:size(temp_h_sciout400,1)
		if strncmp(temp_h_sciout400{n}(6:end),'iPTF',4)==0
			obsdata=dlmread(strcat(char(temp_h_sciout400{n}),num2str(n),'.txt'));
			temp400=textread(strcat(char(temp_h_sci400{n}),num2str(n),'_400.txt'),'%s');
			AirMass=cell2mat(get_fits_keyword(temp400{1},{'AIRMASS'},[]));
			
			EXPTIME=cell2mat(get_fits_keyword(temp400{1},{'EXPTIME'},[]));
			obsdata(:,2)=obsdata(:,2)./EXPTIME;
			
			obsdata=atmospheric_ext(obsdata,AirMass,'mkoextinct.dat','linear'); % atm. ext correction
	
			interpobject=interp1(obsdata(:,1),obsdata(:,2),wl400,'poly1');

	
			flux400std=interpobject.*calfactors';
		end
	end
end


% find sensitivity function
for n=standard_n600;
	obsdata=dlmread(strcat(char(temp_h_sciout600{n}),num2str(n),'.txt'));
	
	temp600=textread(strcat(char(temp_h_sci600{n}),num2str(n),'_600.txt'),'%s');
	AirMass=cell2mat(get_fits_keyword(temp600{1},{'AIRMASS'},[]));
	
	EXPTIME=cell2mat(get_fits_keyword(temp600{1},{'EXPTIME'},[]));
	obsdata(:,2)=obsdata(:,2)./EXPTIME;
			
	obsdata=atmospheric_ext(obsdata,AirMass,'mkoextinct.dat','linear'); % atm. ext correction
	
	interpstandard=interp1(standard.Spec(:,1),standard.Spec(:,2),wl600,'poly1');
	interpobject=interp1(obsdata(:,1),obsdata(:,2),wl600,'poly1');

% % % 	calfactors=smooth(interpstandard,3)./smooth(interpobject,3);  % fit function to exclude tellurics, etc, order 7 polynomial seeems ok, instead of smooth + fit to calfactors
% % % 	calfit=fit(wl600',calfactors,'poly3','normalize','on');
% % % 	calfactors=calfit(wl600);
	
	%exclude atm. abs lines, etc from sens. func. fit
	stdfit=fit(wl600',smooth(interpstandard,3),'poly7','normalize','on');
	objfit=fit(wl600',smooth(interpobject,3),'poly7','normalize','on');
	
	%exclude
		interp_s_hi=interpstandard'<1.2.*stdfit(wl600);
		interp_s_lo=interpstandard'>0.85.*stdfit(wl600);
		interp_s_i=boolean(interp_s_hi.*interp_s_lo);
	
		interp_o_hi=interpobject'<1.2.*objfit(wl600);
		interp_o_lo=interpobject'>0.85.*objfit(wl600);
		interp_o_i=boolean(interp_o_hi.*interp_o_lo);
	
		%interp_i=boolean(interp_o_i.*interp_s_i);
		
		% grow rejection by -10+20 units
		interp_o_i_temp=interp_o_i;
		for i=1:length(interp_o_i)
			if interp_o_i(i)==0 & i > 50 & i < length(interp_o_i)-100
				interp_o_i_temp(i-50:1:i+100)=0;
			end
		end
		
		% grow rejection by -10+20 units
		interp_s_i_temp=interp_s_i;
		for i=1:length(interp_s_i)
			if interp_s_i(i)==0 & i > 50 & i < length(interp_s_i)-100
				interp_s_i_temp(i-50:1:i+100)=0;
			end
		end
	
		interp_o_i=interp_o_i_temp;
		interp_s_i=interp_s_i_temp;
		
		stdfit=fit(wl600(interp_s_i)',smooth(interpstandard(interp_s_i),1),'poly7','normalize','on');
		
		% divide spectrum in 3 parts and fit polynomials

		wl_f=wl600(interp_o_i)';
		interpobject_f=interpobject(interp_o_i);

for m=1:999
		figure(666)
		clf
		hold on
		plot(wl600(interp_o_i),interpobject(interp_o_i),'ko','markersize',2,'markerfacecolor',rgb('SkyBlue'),'color',rgb('SkyBlue'))
		plot(wl600(~interp_o_i),interpobject(~interp_o_i),'ko','markersize',2,'markerfacecolor',rgb('Salmon'),'color',rgb('Salmon'))
		plot(wl600,objfit(wl600),'--','color',rgb('Silver'),'linewidth',1.5)
		title('Sensitiviity function fits, and excluded regions, Observed Spectr. standard')
		box on
		
		x = inputdlg('Enter order of 3 polynomials to fit, e.g. 3 3 3:',...
		'Sample', [1 50]);
		orders = str2num(x{:});	
		
		'Select cutoff points for 3 polynomials to fit'
		
		[input_m,y] = ginput(2);
		input1=input_m(1);
		input2=input_m(2);
		
		wl_f1=wl_f(1:1000+round(mean(wl_f(wl_f<=input1.*1.01 & wl_f >=input1.*0.99))));
		wl_f2=wl_f(-1000+round(mean(wl_f(wl_f<=input1.*1.01 & wl_f >=input1.*0.99))):1000+round(mean(wl_f(wl_f<=input2.*1.01 & wl_f >=input2.*0.99))));
		wl_f3=wl_f(-1000+round(mean(wl_f(wl_f<=input2.*1.01 & wl_f >=input2.*0.99))):end);

		o_f1=smooth(interpobject_f(1:1000+round(mean(wl_f(wl_f<=input1.*1.01 & wl_f >=input1.*0.99)))),1);		
		o_f2=smooth(interpobject_f(-1000+round(mean(wl_f(wl_f<=input1.*1.01 & wl_f >=input1.*0.99))):1000+round(mean(wl_f(wl_f<=input2.*1.01 & wl_f >=input2.*0.99)))),1);
		o_f3=smooth(interpobject_f(-1000+round(mean(wl_f(wl_f<=input2.*1.01 & wl_f >=input2.*0.99))):end),1);
		
		objfit1=fit(wl_f1,o_f1,strcat('poly',num2str(orders(1))),'normalize','on');	
		objfit2=fit(wl_f2,o_f2,strcat('poly',num2str(orders(2))),'normalize','on');	
		objfit3=fit(wl_f3,o_f3,strcat('poly',num2str(orders(3))),'normalize','on');	

		objfit_p1=objfit1(wl600(1):0.2:round(input1));
		objfit_p2=objfit2(round(input1)+0.2:0.2:round(input2));
		objfit_p3=objfit3(round(input2)+0.2:0.2:wl600(end));
		
		stdfit_p=smooth(stdfit(wl600),500);
		objfit_p=smooth([objfit_p1;objfit_p2;objfit_p3],500);

		calfactors=stdfit_p./objfit_p;
	
		figure(666)
		hold on
		plot(wl600,objfit_p,'-','color',rgb('DarkOrchid'),'linewidth',1.5)
		title('Sensitiviity function fits, and excluded regions, Observed Spectr. standard')
		box on
	
		m=m+1;	
		%m = input('If fit is satisfactory type 0 and press enter. Otherwise press enter to redo fit.');
		x = inputdlg('If fit is satisfactory type 0 and press ok. Otherwise press ok to redo fit.',...
		'Sample', [1 50]);
		m = str2num(x{:});	
		if m==0
			break
		end
end

	figure(66)
	subplot(2,3,4)
	hold on
	plot(wl600(interp_s_i),interpstandard(interp_s_i),'ko','markersize',2,'markerfacecolor',rgb('SkyBlue'),'color',rgb('SkyBlue'))
	plot(wl600(~interp_s_i),interpstandard(~interp_s_i),'ko','markersize',2,'markerfacecolor',rgb('Salmon'),'color',rgb('Salmon'))
	plot(wl600,stdfit_p,'-','color',rgb('DarkOrchid'),'linewidth',1.5)
	plot(wl600,stdfit(wl600),'--','color',rgb('Silver'),'linewidth',1.5)

	title('Sensitiviity function fits, and excluded regions, Database Spectr. standard')
	box on

	subplot(2,3,5)
	hold on
	plot(wl600(interp_o_i),interpobject(interp_o_i),'ko','markersize',2,'markerfacecolor',rgb('SkyBlue'),'color',rgb('SkyBlue'))
	plot(wl600(~interp_o_i),interpobject(~interp_o_i),'ko','markersize',2,'markerfacecolor',rgb('Salmon'),'color',rgb('Salmon'))
	plot(wl600,objfit_p,'-','color',rgb('DarkOrchid'),'linewidth',1.5)
	plot(wl600,objfit(wl600),'--','color',rgb('Silver'),'linewidth',1.5)
	title('Sensitiviity function fits, and excluded regions, Observed Spectr. standard')
	box on
	
	subplot(2,3,6)
	plot(wl600,calfactors)
	title('Sensitivity function')
	box on

	% make interactive function to manually choose regions to exclude from standard star for sens. function fit.
	
	
	% apply to science frame(s)
	for n=1:size(temp_h_sciout600,1)
		if strncmp(temp_h_sciout600{n}(6:end),'iPTF',4)==1
			obsdata=dlmread(strcat(char(temp_h_sciout600{n}),num2str(n),'.txt'));
			
			temp600=textread(strcat(char(temp_h_sci600{n}),num2str(n),'_600.txt'),'%s');
			AirMass=cell2mat(get_fits_keyword(temp600{1},{'AIRMASS'},[]));
			
			EXPTIME=cell2mat(get_fits_keyword(temp600{1},{'EXPTIME'},[]));
			obsdata(:,2)=obsdata(:,2)./EXPTIME;
			
			obsdata=atmospheric_ext(obsdata,AirMass,'mkoextinct.dat','linear'); % atm. ext correction
	
			interpobject=interp1(obsdata(:,1),obsdata(:,2),wl600,'poly1');
			
	figure(200+n)

	flux600=interpobject.*calfactors';
			plot(wl600,flux600,'k-','linewidth',1)
			xlabel('Wavelength [Angstr]')
			ylabel('Flux (erg/cm2/s)')
			title('Flux calibrated spectrum')
			box on
		end	
	end
	
		% apply to std frame(s)
	for n=1:size(temp_h_sciout600,1)
		if strncmp(temp_h_sciout600{n}(6:end),'iPTF',4)==0
			obsdata=dlmread(strcat(char(temp_h_sciout600{n}),num2str(n),'.txt'));
			
			temp600=textread(strcat(char(temp_h_sci600{n}),num2str(n),'_600.txt'),'%s');
			AirMass=cell2mat(get_fits_keyword(temp600{1},{'AIRMASS'},[]));
			
			EXPTIME=cell2mat(get_fits_keyword(temp600{1},{'EXPTIME'},[]));
			obsdata(:,2)=obsdata(:,2)./EXPTIME;
			
			obsdata=atmospheric_ext(obsdata,AirMass,'mkoextinct.dat','linear'); % atm. ext correction
	
			interpobject=interp1(obsdata(:,1),obsdata(:,2),wl600,'poly1');
			

			flux600std=interpobject.*calfactors';
	
		end	
	end
end

			figure(999)
			hold on
			plot(wl600,smooth(flux600std,3),'k-','linewidth',1.5,'color',rgb('SteelBlue'))
			plot(wl400,smooth(flux400std,3),'k-','linewidth',1.5,'color',rgb('LightCoral'))
			plot(standard.Spec(:,1),standard.Spec(:,2),'r-','linewidth',1,'color',rgb('Silver'))


			xlabel('Wavelength [Angstr]')
			ylabel('Flux (erg/cm2/s)')
			title('Flux calibrated spectrum')
			box on
			
			

flux_600_temp=smooth(flux600',20)./smooth(flux600,100);
flux600_t=flux600(flux_600_temp<1.2 & flux_600_temp > 0.8);
wl600_t=wl600(flux_600_temp<1.2 & flux_600_temp > 0.8);


flux_400_temp=smooth(flux400',20)./smooth(flux400,100);
flux400_t=flux400(flux_400_temp<1.2 & flux_400_temp > 0.8);
wl400_t=wl400(flux_400_temp<1.2 & flux_400_temp > 0.8);

			figure(1000)
			hold on
			plot(wl600_t,smooth(flux600_t,5),'k-','linewidth',1,'color',rgb('SteelBlue'))
			plot(wl400_t,smooth(flux400_t,5),'k-','linewidth',1,'color',rgb('LightCoral'))


			xlabel('Wavelength [Angstr]')
			ylabel('Flux (erg/cm2/s)')
			title('Flux calibrated spectrum')
			box on
			
			keck=dlmread('keck1.txt');
			plot(keck(:,1),smooth(keck(:,2),3),'k-','linewidth',1,'color',rgb('Silver'))
			
			
			
			
			
			
			
			
			
			
			
			
			
			%%% make average airmass calculation, or read avg airmass from header (start+end + average, or something more correct?)
			%%% make manual selection of areas to exclude from fits
