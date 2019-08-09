function [meanGND,medianGND,stdGND,Vf]= HEDM_GNDDensity(Temperature,Cycle,PixelSpacing,ConfidenceThreshold)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function Definition Line: [meanGND,medianGND,stdGND,Vf]= HEDM_GNDDensity(Temperature,Cycle,PixelSpacing,ConfidenceThreshold)
%  Requires MTEX Version 5.2.2 beta
%  Requires all data files in MATLAB path
%
%  Inputs:
%  1. Temperature (Degrees Celsius) (500, 700, 935, 815, 965, 930, 990, 1000)
%  2. Cycle (1, 2, 3, 4)
%  3. PixelSpacing (microns)
%  4. ConfidenceThreshold (0 <= ConfidenceThreshold <= 1)
%
%  Outputs:
%  1. meanGND: Mean geometrically necessary dislocation density excluding
%  grain boundaries and off sample points (m^(-2))
%  2. medianGND: Median geometrically necessary dislocation density excluding
%  grain boundaries and off sample points (m^(-2))
%  3. stdGND: Standard deviation of geometrically necessary dislocation density excluding
%  grain boundaries and off sample points (m^(-2))
%  4. Vf: Volume fraction of precipitated particles
%
%  Function Description: Estimates the minimum Geometrically Necessary Dislocation
%  Density for HEDM orientation data collected during a cyclical heat
%  treatment experiment at Argonne National Laboratory by Dr. Ashwin Shahani and his graduate students
%  from the University of Michigan 
%  
%  Acknowledgments: 
%  1. MTEX Developers for GND Density codes: http://mtex-toolbox.github.io/team.html
%  2. Procedure: Pantleon W. (2008). Resolving the geometrically necessary
%  dislocation content by conventional electron backscattering diffraction.
%  Scripta Materialia, 58(11), 994-997
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load Data
if (Temperature==500 && Cycle==1);
    fn = 'dummy_2_rt_500_2_nf_bcc_0.25degree_q11_z0_150x150_0.007_shift_0.0_0.0_0.0.h5';
    fn2='dummy_2_rt_500_2_nf_fcc_0.25degree_q11_z0_150x150_0.007_shift_0.0_0.0_0.0.h5';
    ori_test = double(h5read(fn,'/slices/z0/EulerAngles'));
    confi_fcc=h5read(fn2,'/slices/z0/Confidence');
    confi_test = h5read(fn,'/slices/z0/Confidence');
    
elseif (Temperature==700 && Cycle==1);
    fn = 'dummy_2_rt_700_nf_bcc_0.25degree_q11_z0_150x150_0.007_shift_0.0_0.0_0.0.h5';
    fn2='dummy_2_rt_700_nf_fcc_0.25degree_q11_z0_150x150_0.007_shift_0.0_0.0_0.0.h5';
    ori_test = double(h5read(fn,'/slices/z0/EulerAngles'));
    confi_fcc=h5read(fn2,'/slices/z0/Confidence');
    confi_test = h5read(fn,'/slices/z0/Confidence');
    
elseif (Temperature==935 && Cycle==1)
    fn = 'dummy_2_rt_935_nf_bcc_0.25degree_q11_z0_150x150_0.007_shift_0.0_0.0_0.0.h5';
    fn2='dummy_2_rt_935_nf_fcc_0.25degree_q11_z0_150x150_0.007_shift_0.0_0.0_0.0.h5';
    ori_test = double(h5read(fn,'/slices/z0/EulerAngles'));
    confi_fcc=h5read(fn2,'/slices/z0/Confidence');
    confi_test = h5read(fn,'/slices/z0/Confidence');
    
elseif (Temperature==500 && Cycle==2);
    fn = 'dummy_2_rt_500_nf_restart_bcc_0.25degree_q11_z1_150x150_0.007_shift_0.0_0.0_0.0.h5';
    fn2='dummy_2_rt_500_nf_restart_fcc_0.25degree_q11_z1_150x150_0.007_shift_0.0_0.0_0.0.h5';
    ori_test = double(h5read(fn,'/slices/z1/EulerAngles'));
    confi_fcc=h5read(fn2,'/slices/z1/Confidence');
    confi_test = h5read(fn,'/slices/z1/Confidence');
    
elseif (Temperature==700 && Cycle==2);
    fn = 'dummy_2_rt_700_2_nf_restart_2_bcc_0.25degree_q11_z0_150x150_0.007_shift_0.0_0.0_0.0.h5';
    fn2='dummy_2_rt_700_2_nf_restart_2_fcc_0.25degree_q11_z0_150x150_0.007_shift_0.0_0.0_0.0.h5';
    ori_test = double(h5read(fn,'/slices/z0/EulerAngles'));
    confi_fcc=h5read(fn2,'/slices/z0/Confidence');
    confi_test = h5read(fn,'/slices/z0/Confidence');
    
elseif (Temperature==815 && Cycle==2);
    fn = 'dummy_2_rt_815_2_nf_copperBCC_q11_rot720_26layers_150x150_0.007_shift_0.0_0.0_0.0.h5';
    fn2='dummy_2_rt_815_2_nf_copperFCC_q11_rot720_26layers_150x150_0.007_shift_0.0_0.0_0.0.h5';
    ori_test = double(h5read(fn,'/slices/z19/EulerAngles'));
    confi_fcc=h5read(fn2,'/slices/z19/Confidence');
    confi_test = h5read(fn,'/slices/z19/Confidence');
    
elseif (Temperature==965 && Cycle==2);
    fn = 'dummy_2_single_crystal_furnace_nf_bcc_0.25degree_q11_z1_150x150_0.007_shift_0.0_0.0_0.0.h5';
    fn2='dummy_2_single_crystal_furnace_nf_fcc_0.25degree_q11_z1_150x150_0.007_shift_0.0_0.0_0.0.h5';
    ori_test = double(h5read(fn,'/slices/z1/EulerAngles'));
    confi_fcc=h5read(fn2,'/slices/z1/Confidence');
    confi_test = h5read(fn,'/slices/z1/Confidence');
    
elseif (Temperature==930 && Cycle==3);
    fn = 'dummy_2_cycle3_rt_930_furnace_nf_bcc_0.25degree_q11_z0_150x150_0.007_shift_0.0_0.0_0.0.h5';
    fn2='dummy_2_cycle3_rt_930_furnace_nf_fcc_0.25degree_q11_z0_150x150_0.007_shift_0.0_0.0_0.0.h5';
    ori_test = double(h5read(fn,'/slices/z0/EulerAngles'));
    confi_fcc=h5read(fn2,'/slices/z0/Confidence');
    confi_test = h5read(fn,'/slices/z0/Confidence');
    
elseif (Temperature==990 && Cycle==3);
    fn = 'dummy_2_cycle3_rt_990_furnace_nf_bcc_0.25degree_q11_z1_150x150_0.007_shift_0.0_0.0_0.0.h5';
    fn2='dummy_2_cycle3_rt_990_furnace_nf_fcc_0.25degree_q11_z1_150x150_0.007_shift_0.0_0.0_0.0.h5';
    ori_test = double(h5read(fn,'/slices/z1/EulerAngles'));
    confi_fcc=h5read(fn2,'/slices/z1/Confidence');
    confi_test = h5read(fn,'/slices/z1/Confidence');
    
elseif (Temperature==1000 && Cycle==4);
    fn = 'dummy_2_cycle4_rt_1000_furnace_nf_bcc_0.25degree_q11_z1_150x150_0.007_shift_0.0_0.0_0.0.h5';
    fn2='dummy_2_cycle4_rt_1000_furnace_nf_fcc_0.25degree_q11_z1_150x150_0.007_shift_0.0_0.0_0.0.h5';
    ori_test = double(h5read(fn,'/slices/z1/EulerAngles'));
    confi_fcc=h5read(fn2,'/slices/z1/Confidence');
    confi_test = h5read(fn,'/slices/z1/Confidence');
    
else
    error('Please select a valid [Temperature Cycle] combination: [500 1], [700 1], [935 1], [500 2], [700 2], [815 2], [965 2], [930 3], [990 3], [1000 4]')
end
Temperature=num2str(Temperature);
Cycle=num2str(Cycle);

%Symmetries and Crystalographic Information
CS=crystalSymmetry('432',[2.947 2.947 2.947], [90 90 90]*degree);
SS=specimenSymmetry('1');
a=norm(CS.aAxis);

%Create Orientation Objects
[l,m,n]=size(ori_test);
ori=orientation.byEuler((reshape(ori_test(1,:,:),m,n)),(reshape(ori_test(2,:,:),m,n)),(reshape(ori_test(3,:,:),m,n)),CS,SS);

%Define DislocationSystem and Properties
dS=dislocationSystem.bcc(CS);
nu=0.3;
dS(dS.isEdge).u=1;
dS(dS.isScrew).u=1-nu;
dSRot=ori*dS;

%Data Parameters

dx=PixelSpacing;
dy=dx;

%Get Curvatures

    %MTEX Retrieved Code (EBSDsquare.m)
    ori_right = ori(:,[2:end end-1]);
    gX = log(ori_right,ori,'left') ./ dx;
    gX(:,end) = - gX(:,end);
    ori_up = ori([2:end end-1],:);
    gY = log(ori_up,ori,'left') ./ dy;
    gY(end,:) = - gY(end,:);
    
    %MTEX Retrieved Code (curvature.m)
    kappa = dyad(gX,tensor([1;0;0])) + ...
        dyad(gY,tensor([0;1;0]));
    kappa{:,3}=NaN;
    kappa = curvatureTensor(kappa,'unit',['1/um']);
   
%Caculate Dislocation Density
[rho,factor]=fitDislocationSystems(kappa,dSRot);
alpha=sum(dSRot.tensor .* rho,2);
alpha.opt.unit='1/um';
gndtotal=reshape(factor*sum(abs(rho .* dSRot.u),2),m,n)';

%Apply Confidence Filter
confi=zeros(m,n);
for idxm=1:m;
    for idxn=1:n;
        if confi_test(idxm,idxn) < ConfidenceThreshold | confi_test(idxm,idxn) > 1;
           confi(idxm,idxn)=0;
        else
            confi(idxm,idxn)=confi_test(idxm,idxn);
        end
        if confi_fcc(idxm,idxn) < ConfidenceThreshold | confi_fcc(idxm,idxn) > 1;
           confi_fcc(idxm,idxn)=0;
        end
        
        %Superimpose FCC Particles on BCC Matrix
        if confi(idxm,idxn)< confi_fcc(idxm,idxn);
           confi(idxm,idxn)=NaN;
        elseif confi(idxm,idxn)~=0;
           confi(idxm,idxn)=1;
        end
    end
end
gndtotal=gndtotal.*confi;

%Plot minimum geometrically necessary dislocation density
contourf(gndtotal,1000,'edgecolor','none')
CLim=[1e12 1e14];
set(gca,'ColorScale','log');
set(gca,'CLim',CLim);
axis square;
colorbar;
title(sprintf('Minimum Geometrically Necessary Dislocation Density (1/m^{2}) Cycle %s %s%cC',Cycle,Temperature,char(176)));

%Statistical Analysis
%Volume Fraction Particles
Vf=sum(isnan(confi)==1,'all')/(m*n);

%Off grain boundary GND Densities
gndtotalNoGB=double((gndtotal>(0) & gndtotal<(max(CLim,[],'all'))));
gndtotalNoGB(gndtotalNoGB==0)=NaN;
gndtotalNoGB=gndtotalNoGB.*gndtotal;
meanGND=mean(gndtotalNoGB,'all','omitnan');
medianGND=median(gndtotalNoGB,'all','omitnan');
stdGND=std(gndtotalNoGB,0,'all','omitnan');
fprintf('Measures of central tendency and variation of GND Density in (1/m^(2))\nMean: %g\nMedian: %g\nStandard Deviation: %g\n', round(meanGND,3,'significant'), round(medianGND,3,'significant'), round(stdGND,3,'significant'))
fprintf('Estimated volume fraction of particles: %g\n',round(Vf,3,'significant'))
