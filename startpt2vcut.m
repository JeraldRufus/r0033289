%% identify range of ref pts from start pt
% input ref pts and start pts
ref_pts= xlsread('midref_and_ref_pts.xlsx',3, 'P:R');
startpt_data=readtable('SL LIGHT GUIDE _210222 _0.5CS_OPT2.txt');
start_pts=startpt_data(:,2:4);
start_pts=start_pts{:,:};   %coverts table to array
start_dirvec=startpt_data(:,5:7);
start_dirvec=start_dirvec{:,:};

%input source pt for initial sorting
source_pt= [-130.52 170.005 2.268];
for i=1:size(start_pts) 
    d(i,:)=[i, sqrt((source_pt(1,1)- start_pts(i,1))^2 + (source_pt(1,2)- start_pts(i,2))^2 + (source_pt(1,3)- start_pts(i,3))^2), start_pts(i,1),start_pts(i,2),start_pts(i,3)];  
  
end
k=1;
for i=1:size(d,1)
   if d(i,2)>15
       d1(k,:)=d(i,:); % d1 is new ray set
       k=k+1;
   end
end   
k=1;p=1; %% corresponds to start pt falling within 0 to 5mm distance of ref pt
for i=1:size(d1)
    for j=1:size(ref_pts)
    pts_dist(j,i)=abs(sqrt((ref_pts(j,1)- d1(i,3))^2 + (ref_pts(j,2)- d1(i,4))^2 + (ref_pts(j,3)- d1(i,5))^2));
    
    end
%     if min(pts_dist(:,i))< 15  %%ask
%         pts_dist(:,i)=pts_dist(:,i);
    idx(k)=find(pts_dist(:,i)==min(pts_dist(:,i)));
    idx2(k)=idx(k)-1; idx3(k)=idx(k)+1;
        if pts_dist(idx2(k))>pts_dist(idx3(k))
            idx_new(k)=idx3(k);
        else
            idx_new(k)=idx2(k);
        end
    vcut(k)=idx(k);vcut2(k)=idx_new(k); ray(k)=d1(i,1); k=k+1;
%     else
%         pts_dist(:,i)=zeros(size(pts_dist(:,i)));
%         d2(p)= d1(i,1);p=p+1;
%     end
    
end
% d2=d2'; %% start pts at >15mm distance from refpts
vcut2ray=[vcut' ray'];
refpts2startpt=[vcut', vcut2', ray'];
xlswrite('startpt2vcut.xlsx', vcut2ray) %%change filename

%% creation of surface normals at each start point
ref_norm= xlsread('midref_and_ref_pts.xlsx',3, 'V:W'); n=size(ref_norm,1);
ref_norm=[ (linspace(1,n,n))', ref_norm];

% for i=1:size(refpts2startpt,1)
%    var1=refpts2startpt(i,1); var2=refpts2startpt(i,2);
%    var1_norm=ref_norm(var1,2:4); var1_ref=ref_pts(var1,:);
%    var2_norm=ref_norm(var1,2:4); var2_ref=ref_pts(var2,:);
%    dir=(var2_ref(1,1)-var1_ref(1,1))^2 + (var2_ref(1,2)- var1_ref(1,2))^2+ (var2_ref(1,3)- var1_ref(1,3))^2;
%    dir_d=[ var2_ref(1,1)-var1_ref(1,1) , var2_ref(1,2)- var1_ref(1,2) , var2_ref(1,3)- var1_ref(1,3)];
%    vec_d= dir_d/sqrt(dir);
%    var3=refpts2startpt(i,3); var3_st=start_pts(var3,:);
%    vec_v=[var3_st(1,1)-var1_ref(1,1) , var3_st(1,2)- var1_ref(1,2) , var3_st(1,3)- var1_ref(1,3)];
%    vec_t=dot(vec_v,vec_d);
%    new_p(i,:)= var1_ref + vec_t * vec_d;
%    x_norm= var1_norm(1,1) + ( new_p(i,1) -var1_ref(1,1)) * (var2_norm(1,1)-var1_norm(1,1) /(var2_ref(1,1)-var1_ref(1,1)));
%    y_norm= var1_norm(1,2) + ( new_p(i,2) -var1_ref(1,2)) * (var2_norm(1,2)-var1_norm(1,2) /(var2_ref(1,2)-var1_ref(1,2)));
%    z_norm= var1_norm(1,3) + ( new_p(i,3) -var1_ref(1,3)) * (var2_norm(1,3)-var1_norm(1,3) /(var2_ref(1,3)-var1_ref(1,3)));
%    norm_vec(i,:)= [x_norm, y_norm, z_norm];
% 
% 
% end

rib_height=0.702; %%user input
lg_radius=3; %%user input
t=rib_height + lg_radius;

for i=1:size(refpts2startpt,1)
   var1=refpts2startpt(i,1); var2=refpts2startpt(i,2);
   var1_norm=ref_norm(var1,2:4); var1_ref=ref_pts(var1,:);
   var2_norm=ref_norm(var1,2:4); var2_ref=ref_pts(var2,:);
   
   dir=(var2_ref(1,1)-var1_ref(1,1))^2 + (var2_ref(1,2)- var1_ref(1,2))^2+ (var2_ref(1,3)- var1_ref(1,3))^2;
   dir_d=[ var2_ref(1,1)-var1_ref(1,1) , var2_ref(1,2)- var1_ref(1,2) , var2_ref(1,3)- var1_ref(1,3)]; %%change
   vec_d= dir_d/sqrt(dir);
   
   var3=refpts2startpt(i,3); var3_st=start_pts(var3,:);
   
   vec_v=[var3_st(1,1)-var1_ref(1,1) , var3_st(1,2)- var1_ref(1,2) , var3_st(1,3)- var1_ref(1,3)]; %%change
   
   v_t=dot(vec_v,vec_d);
   new_p(i,:)= var1_ref + v_t * vec_d;
   
   dref1=sqrt((var1_ref(1,1)- new_p(i,1))^2 + (var1_ref(1,2)- new_p(i,2))^2 + (var1_ref(1,3)- new_p(i,3))^2);  
   dref2=sqrt((var2_ref(1,1)- new_p(i,1))^2 + (var2_ref(1,2)- new_p(i,2))^2 + (var2_ref(1,3)- new_p(i,3))^2); 
   dref=sqrt((var1_ref(1,1)- var2_ref(1,1))^2 + (var1_ref(1,2)- var2_ref(1,2))^2 + (var1_ref(1,3)- var2_ref(1,3))^2); 
   
   x_norm= (var1_norm(1,1)*dref1 + var2_norm(1,1)* dref2)/dref;
   y_norm= (var1_norm(1,2)*dref1 + var2_norm(1,2)* dref2)/dref;
   z_norm= (var1_norm(1,3)*dref1 + var2_norm(1,3)* dref2)/dref;
   
   x_n= x_norm/sqrt((x_norm)^2 + (y_norm)^2+ (z_norm)^2);
   y_n= y_norm/sqrt((x_norm)^2 + (y_norm)^2+ (z_norm)^2);
   z_n= z_norm/sqrt((x_norm)^2 + (y_norm)^2+ (z_norm)^2);
   norm_vec(i,:)= [x_n, y_n, z_n];
   
   center_dist(i,:) = new_p(i,:)+ t *norm_vec(i,:);
   surf_norm(i,:)= [var3, (var3_st - center_dist(i,:))/sqrt((var3_st(1,1)-center_dist(i,1))^2 + (var3_st(1,2)- center_dist(i,2))^2+ (var3_st(1,3)- center_dist(i,3))^2)];
  

end

%% Backward ray tracing from start point to v-cut
ref_ind_lg=1.49; %user input
ref_ind_air=1; %user input
for i=1:size( surf_norm,1)
    a_unit=surf_norm(i,2:4); %unit vector, mag=1
    var4=surf_norm(i,1); refr_rayvec=start_dirvec(var4,:); mag_refr_rayvec=sqrt(refr_rayvec(1,1)^2 + refr_rayvec(1,2)^2 + refr_rayvec(1,3)^2);
    ang_r(i,:)= [var4, acosd((dot(a_unit,refr_rayvec))/( mag_refr_rayvec))];
    r_t=dot(a_unit,refr_rayvec); a_vec=r_t * a_unit;
    b_vec= refr_rayvec - a_vec;
    b_unit= b_vec/ sqrt(b_vec(1,1)^2 + b_vec(1,2)^2 + b_vec(1,3)^2);
    ang_i(i,:)= [var4, asind((ref_ind_air/ref_ind_lg) * sind(ang_r(i,2)))];
    inci_rayvec(i,:)=  [var4, (cosd(ang_i(i,2)) * a_unit +  sind(ang_i(i,2)) * b_unit)];

    
    
end


%% creation of planes at each v-cut
% ref_pts= xlsread('DOE_0402_3_0.5mmcontHeight_213.xlsx',2, 'B:D');
% norm_vec=xlsread('DOE_0402_3_0.5mmcontHeight_213.xlsx',1, 'AZ:BB');
% 
% 
% %% creation on bounding boxes at each v-cut
% midref_pts= xlsread('DOE_0402_3_0.5mmcontHeight_213 - midrefpts.xlsx',2, 'B:D');
% norm_vec=xlsread('DOE_0402_3_0.5mmcontHeight_213 - midrefpts.xlsx',1, 'AZ:BB');
% line_vec_start=[]; line_vec_end=[];
% for i=1:size(midref_pts)
%     line_vec_start(i,:)=midref_pts(i,:)+ 1 * width_vec(i,:);
% 
% end
% for i=1:size(ref_pts)
%     line_vec_end(i,:)=midref_pts(i,:)+ (-1) * width_vec(i,:);
% end
% 
% %% project boxes on the v-cut plane
% 
% 
