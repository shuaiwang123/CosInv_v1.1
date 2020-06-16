% This script is used to find the edges of fault patches, and finally
% return the indexes of the patches who lie in the edge.
% close all

eps=1.e-1;

% x0=fault_model(:,1);y0=fault_model(:,2);z0=fault_model(:,3);
% fault_file=load('../Coupling_Peru/faultmodel/peru.trg');
%
%---fault_file='faultmodel/nias.trg'
%---x0,y0,z0断层片的中心在局部坐标系中的坐标[单位：Km]
x0=fault_file(:,1);y0=fault_file(:,2);z0=fault_file(:,3);
xc=0.;yc=0.;
theta = strike_average;
%---rotate2d_center_matrix    /Code/Fault Related
%   以xc,yc为中心,对x0,y0点逆时针旋转一个theta进行坐标转换,并输出转换后的坐标
[Xtmp,Ytmp]=rotate2d_center_matrix(x0,y0,xc,yc,theta);
z=z0;

npts=numel(Xtmp);

% figure(1);
% plot3(x0,y0,-z0,'.b');
% 
% dx=0;dy=dx;
% 
% hold on
% for ii=1:npts,
%     text(x0(ii)+dx,y0(ii)+dy,num2str(ii));
% end
% view(2)


x_min=min(Xtmp);x_max=max(Xtmp);
y_min=min(Ytmp);y_max=max(Ytmp);

dy1=sort(abs(diff(Ytmp)));
dx1=sort(abs(diff(Xtmp)));

ind_y1=find(dy1>eps);dy1=dy1(ind_y1);dy=dy1(1);
ind_x1=find(dx1>eps);dx1=dx1(ind_x1);dx=dx1(1);

n_int_y=length(y_min:dy:y_max);

% plot3(Xtmp,Ytmp,-z,'.b');view(2);
% hold on;

yi=y_min-dy/2;
yf=yi+dy;
right_point=zeros(3,n_int_y);
left_point=zeros(3,n_int_y);

index_lr=zeros(2*n_int_y,1);

for ii=1:n_int_y,
    ic=find(Ytmp>=yi & Ytmp<=yf);
    x_tmp=Xtmp(ic);
    
    if ~isempty(x_tmp),
        [x_tmp_sorted,ind_x_tmp_sorted]=sort(x_tmp,'ascend');
        ind_left_pt=find(Xtmp==x_tmp_sorted(1));
        ind_right_pt=find(Xtmp==x_tmp_sorted(end));
        
        ind_left_pt=intersect(ic,ind_left_pt);
        ind_right_pt=intersect(ic,ind_right_pt);
        
        left_point(1,ii)=Xtmp(ind_left_pt);
        left_point(2,ii)=Ytmp(ind_left_pt);
        left_point(3,ii)=z(ind_left_pt);
        
        right_point(1,ii)=Xtmp(ind_right_pt);
        right_point(2,ii)=Ytmp(ind_right_pt);
        right_point(3,ii)=z(ind_right_pt);
        
        index_lr(2*ii-1)=ind_left_pt;
        index_lr(2*ii)=ind_right_pt;
        
%         plot3(Xtmp(ind_left_pt),Ytmp(ind_left_pt),-z(ind_left_pt),'.m');
%         plot3(Xtmp(ind_right_pt),Ytmp(ind_right_pt),-z(ind_right_pt),'.m');
%         view(2);
        
    end
    yi=yf;
    yf=yi+dy;
end

index_lr=index_lr(find(index_lr~=0));
index_lr=unique(index_lr);

% close all
% figure(1);
% plot3(Xtmp,Ytmp,-z,'.b');
% hold on;
% plot3(left_point(1,:),left_point(2,:),-left_point(3,:),'.m');
% plot3(right_point(1,:),right_point(2,:),-right_point(3,:),'.m');
% view(2);
% keyboard

xi=x_min-dx/2;
xf=xi+dx;

n_int_x=length(x_min:dx:x_max);

upper_point=zeros(3,n_int_x);
lower_point=zeros(3,n_int_x);

index_ud=zeros(2*n_int_x,1);

for ii=1:n_int_x,
    ic=find(Xtmp>=xi & Xtmp<=xf);
    y_tmp=Ytmp(ic);
    
    if ~isempty(y_tmp),
        [y_tmp_sorted,ind_y_tmp_sorted]=sort(y_tmp,'ascend');
        ind_lower_pt=find(Ytmp==y_tmp_sorted(1));
        ind_upper_pt=find(Ytmp==y_tmp_sorted(end));
        
        ind_lower_pt=intersect(ic,ind_lower_pt);
        ind_upper_pt=intersect(ic,ind_upper_pt);
        
        lower_point(1,ii)=Xtmp(ind_lower_pt);
        lower_point(2,ii)=Ytmp(ind_lower_pt);
        lower_point(3,ii)=z(ind_lower_pt);
        
        upper_point(1,ii)=Xtmp(ind_upper_pt);
        upper_point(2,ii)=Ytmp(ind_upper_pt);
        upper_point(3,ii)=z(ind_upper_pt);
                
        index_ud(2*ii-1)=ind_lower_pt;
        index_ud(2*ii)=ind_upper_pt;
%                 keyboard
%         plot3(Xtmp(ind_lower_pt),Ytmp(ind_lower_pt),-z(ind_lower_pt),'.r');
%         plot3(Xtmp(ind_upper_pt),Ytmp(ind_upper_pt),-z(ind_upper_pt),'.r');
    end
    %     keyboard
    xi=xf;
    xf=xi+dx;
end
index_ud=index_ud(find(index_ud~=0));
index_ud=unique(index_ud);

% plot3(lower_point(1,:),lower_point(2,:),-lower_point(3,:),'.r');
% plot3(upper_point(1,:),upper_point(2,:),-upper_point(3,:),'.r');
% view(2);
% keyboard

ind_patches_edge=[index_lr;index_ud];
ind_patches_edge=unique(ind_patches_edge);
ind_patches_edge=ind_patches_edge';

% plot3(x0(ind_patches_edge),y0(ind_patches_edge),-z0(ind_patches_edge),'sr');view(2);

% save ind_patches_edge.mat ind_patches_edge;
% keyboard
% figure(gcf+1)
% plot_labeled_points