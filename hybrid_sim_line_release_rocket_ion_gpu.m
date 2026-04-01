%% calculate the total electron contact of ionosphere model
%%set normal var
clc;
clear;
tic
%iri电子密度模型路径
filename = 'F:\Simulation_rocket_release_ion\hybrid_simulation\irimodel\iri2014_142_ut14.txt';lt=10;
% 打开文件
kk=1.380649 * 10^-23;  %玻尔兹曼常数J/k
gg=9.81 ;%m/s2 1
e0=1.6*10^-19;%电子带电量e库伦 C1.602176634
NA=6.02*10^23;%阿伏伽德罗常数
mh2o=18 ;%H2O分子质量
mh2=2   ;%H2
k1=3.2*10^-9 ; %水反应e的损失系数，单位cm2/s
k12=6.5*10^-7 ; %水反义步骤2，需要再乘以（300/Te）^0.5
k2=2.0*10^-9 ; %h2反应e的损失系数，单位cm2/s
pi=3.1415926536;
%RL-10C-1发动机喷射质量速度为23.078kg/s；反应混合质量比为5.5:1
%水质量占比为198/208，氢气占比质量为10/208
mv0=23.078;%单位kg/s
NvH2o=mv0*10^3 * 198/208 / 18;%水释放速度，单位mol/s
NvH2=mv0*10^3 * 10/208 /2;%氢气的释放速度，单位mol/s
B_abs=37513.1*10^-9;%磁场强度
% B_inc=-53.6926/180*pi;%磁场倾角
% B_sit=-14.91022/180*pi;%磁场转角
load('B_d.mat');
load('B_i.mat');
B_ii=B_i*(pi./-180);
B_di=B_d*(pi./180);
B_inc=gpuArray(B_ii);
B_sit=gpuArray(B_di);
%简化计算
B_i_cos=cos(B_inc);B_i_sin=sin(B_inc);
B_d_cos=cos(B_sit);B_d_sin=sin(B_sit);
B1_a=cos(B_inc).*sin(B_inc).*cos(B_sit);
B2_a=cos(B_inc).*sin(B_inc).*sin(B_sit);
B3_a=sin(B_inc).^2;
clear B_i B_ii B_d B_di B_inc B_sit;
%f.e20.FXSD.f19_f19.001.cam.h1.2014-05-22-00000.nc
%% 读取IRI2016模型数据
fid = fopen(filename, 'r');
if fid == -1
    error('无法打开文件 %s', filename);
end
e_alt = []; % 高度
e_density = []; % 电子密度
T_e=[];%电子温度
o_ion = []; % O+ 百分比*10
T_air=[];%大气温度
T_ion=[];%离子温度
while ~feof(fid)
    line = fgetl(fid); % 读取一行
    if ischar(line) % 检查是否为有效行
        % 用空格分割字符串，并将结果转换为数字
        values = str2num(line); %#ok<ST2NM>
        if ~isempty(values) % 确保有有效数据
            e_alt(end + 1) = values(1); % 第一列是高度
            e_density(end + 1) = values(2); % 第二列是电子密度njhm,.
            T_e(end+1) = values(6); %电子温度
            o_ion(end+1) = values(7); % 氧离子百分比*10
            T_air(end+1) = values(4); % 大气温度
            T_ion(end+1) = values(5); % 离子温度
        end
    end
end
fclose(fid);% 关闭文件

%%MSIS00大气模型参数获取
year=2014;doy=142;sec=14*3600;alt=300;lat=30;lon=-70;lst=12;%alt=300km高度
f107_iri_nc=ncread("f.e20.FXSD.f19_f19.001.cam.h1.2014-05-22-00000.nc",'f107');
f107a_iri_nc=ncread("f.e20.FXSD.f19_f19.001.cam.h1.2014-05-22-00000.nc",'f107a');
ap_iri_nc=ncread("f.e20.FXSD.f19_f19.001.cam.h1.2014-05-22-00000.nc",'ap');%8->5
ap_a=mean(ap_iri_nc);
% density=nrlmsise00_own(year,doy,sec,alt,lat,lon,lst,f107_iri_nc(1),ap_iri_nc(5));
% den_o=density(1).d(2);% O 原子密度
% den_o2=density(1).d(4);% O2分子密度
% den_N2=density(1).d(3);% N2分子密度
% den_H=density(1).d(7);% H原子
% den_N=density(1).d(8);% N原子

Release_Va=101.81*10^3/23.078;%释放物速度，单位m/s
% T0 = T_4D(117,64,32,1);%200km,Z3->1-145
%设置网格大小
%数组e获取电子背景密度;数组o获取O+密度
A_e=gpuArray(zeros(300,110,200));
A_o_ion=A_e;
A_e_T=A_e;
A_h2o_ion=A_e;
A_wa=A_e;%碰撞频率
for i=1:1:200
    A_e(:,:,i)=A_e(:,:,i)+e_density(i+50);
    A_e_T(:,:,i)=A_e_T(:,:,i)+T_e(i+50);
    if o_ion(i+50)==0%氧原子离子含量=0导致出错，设定一个小值优化他且不会影响模型的结果
        o_ion(i+50)=1;%这里是千分比,设置一个低值让他不为零
    end
    A_o_ion(:,:,i)=A_o_ion(:,:,i)+o_ion(i+50)*e_density(i+50)*0.001;%0.001是对o成分千分比修正
end
%% %数组水含量x,y,z%设置范围大小
Array_water=gpuArray(zeros(300,110,200));
%Array-water-velocity-z，水的速度分量z
A_w_v_z=Array_water;%速度分量，z

%扩散系数
D_h2o=Array_water+1;
D_h2o_ion=Array_water+1;
%获取不同高度上各项粒子密度参数
den_o_z=1;den_o2_z=1;den_N2_z=1;
T_a=gpuArray(zeros(300,110,200));
T_i=T_a;A_ma=T_a;A_q=T_a;H_h2o_ion=T_a;
for i=1:1:200 %高度范围
    den_i=nrlmsise00_own(year,doy,sec,100+i*2,lat,lon,lst,f107_iri_nc(1),ap_iri_nc(5));
    den_o_z=den_i(1).d(2);
    den_o2_z=den_i(1).d(4);
    den_N2_z=den_i(1).d(3);
%     if i<150
        D_h2o(:,:,i)=D_h2o(:,:,i).*Dh2o(den_o_z,den_N2_z,den_o2_z,T_air(50+i));%H2o的扩散系数
%     else
%         D_h2o(:,:,i)=D_h2o(:,:,i).*Dh2o(den_o_z,den_N2_z,den_o2_z,T_air(150));
%     end
%     if i<150
        D_h2o_ion(:,:,i)=D_h2o_ion(:,:,i).*Dh2o_i(den_o_z,den_N2_z,den_o2_z,T_ion(50+i)); %h2o+扩散系数
%     else
%         D_h2o_ion(:,:,i)=D_h2o_ion(:,:,i).*Dh2o_i(den_o_z,den_N2_z,den_o2_z,T_ion(150)); %h2o+扩散系数
%     end
    n_a=sum(den_i(1).d([1,2,3,4,5,7,8]));%分子数n_a,cm-3
    m_a=den_i(1).d(6)/n_a;%碰撞粒子平均质量
    A_ma(:,:,i)=A_ma(:,:,i)+m_a;%平均分子质量
    delta_0=5*10^-15;%碰撞截面,cm-2
    A_wa(:,:,i)=n_a*delta_0*(8*kk*T_e(i+50)*(m_a*NA+mh2o)/m_a/mh2o/pi)^0.5;%碰撞频率
    T_a(:,:,i)=T_a(:,:,i)+T_air(i+50);% air温度
    T_i(:,:,i)=T_i(:,:,i)+T_ion(i+50);% ion温度
    H_h2o_ion(:,:,i)=H_h2o_ion(:,:,i)+(T_air(i+50)+T_ion(i+50))*kk*NA/mh2o/gg*1000;%h2o+的大气标高,单位m
    solar_zenith_angle=solarze(lt,142,lat);%太阳天顶角
    e_q_h=epro(6000,120,m_a,100+i*2,T_air(i+50),solar_zenith_angle);%光电离生成率
    if e_q_h<0
        e_q_h=0;
    end
    A_q(:,:,i)=A_q(:,:,i)+e_q_h;%生成率
end

%%创建一个cell用以存储迭代计算结果,t
Cell_water=cell(1,100);Cell_h2o_ion=cell(1,100);Cell_B1=cell(1,100);Cell_B2=cell(1,100);Cell_B3=cell(1,100);Cell_B_D=cell(1,100);
Cell_e=cell(1,100);Cell_o=cell(1,100);%cell_h2o_l={};cell_h2o_r={};cell_h2o_v={};cell_h2o_d={};cell_xpp={};cell_ypp={};cell_zpp={};
[lenx,leny,lenz]=size(Array_water);

%步长设置dt/dx << 1
dx=10;%单位1km
dy=10;
dz=2 ;
dt=0.01;%时间步长 秒,/s
%% %%微分方程计算n
t0=0;kkk=0;
tend=dt*3000*20; %20mins=60000计算120s模型需要时间400秒 1:3计算耗时
while t0<tend
    if t0==0 %初始化偏导参数
        Array_water_zp=Array_water;%水的参数
        Array_water_xpp=Array_water;
        Array_water_ypp=Array_water;
        Array_water_zpp=Array_water;
        A_w_v_z=Array_water;
        A_h2o_i_xpp=A_h2o_ion;%h2o+参数
        A_h2o_i_ypp=A_h2o_ion;
        A_h2o_i_zpp=A_h2o_ion;
        A_h2o_i_xp=A_h2o_ion;
        A_h2o_i_yp=A_h2o_ion;
        A_h2o_i_zp=A_h2o_ion;
        T_i_zp=T_i;%温度的垂直变化一阶微分
        T_i_zp(:,:,2:lenz)=(T_i(:,:,2:lenz)-T_i(:,:,1:lenz-1))/dz*10^-5;
        T_i_zpp=T_i;%温度垂直变化的二阶微分
        T_i_zpp(:,:,2:lenz-1)=(T_i(:,:,1:lenz-2)-2*T_i(:,:,2:lenz-1) ...
            +T_i(:,:,3:lenz))/(dz*dz)*10^-10;
        H_h2o_i_zp=H_h2o_ion;%大气标高的z轴一阶微分
        H_h2o_i_zp(:,:,2:lenz)=(H_h2o_ion(:,:,2:lenz)-H_h2o_ion(:,:,1:lenz-1))/dz*10^-3;
        D_h2o_i_zp=D_h2o_ion;%扩散系数的z轴一阶微分
        D_h2o_i_zp(:,:,2:lenz)=(D_h2o_ion(:,:,2:lenz)-D_h2o_ion(:,:,1:lenz-1))/dz*10^-5;
        %初始化参数o的损失系数
        B_o_los=A_q./A_o_ion;
        B_o_los(isinf(B_o_los))=A_q(isinf(B_o_los));
        B_o_los(isnan(B_o_los))=A_q(isnan(B_o_los));
        D_H=D_h2o_ion./H_h2o_ion;
        B3_a_1=D_h2o_i_zp + D_h2o_ion./T_i.*T_i_zp+D_H;
    end
    %% 计算n对xyz一阶偏导
    %一阶偏导
    Array_water_zp(:,:,2:lenz)=( ...
        Array_water(:,:,2:lenz)-Array_water(:,:,1:lenz-1))/dz;
    A_h2o_i_xp(2:lenx,:,:)=( ...
        A_h2o_ion(2:lenx,:,:)-A_h2o_ion(1:lenx-1,:,:))/dx ;
    A_h2o_i_yp(:,2:leny,:)=( ...
        A_h2o_ion(:,2:leny,:)-A_h2o_ion(:,1:leny-1,:))/dy;
    A_h2o_i_zp(:,:,2:lenz)=( ...
        A_h2o_ion(:,:,2:lenz)-A_h2o_ion(:,:,1:lenz-1))/dz;
    %% 分别计算xyz二阶偏微分
    %x的二阶偏导
    Array_water_xpp(2:lenx-1,:,:)=( ...
        Array_water(3:lenx,:,:)-2*Array_water(2:lenx-1,:,:)+ ...
        Array_water(1:lenx-2,:,:))/(dx*dx);
    A_h2o_i_xpp(2:lenx-1,:,:)=( ...
        A_h2o_ion(3:lenx,:,:)-2*A_h2o_ion(2:lenx-1,:,:)+ ...
        A_h2o_ion(1:lenx-2,:,:))/(dx*dx);
    %y的二阶偏导
    Array_water_ypp(:,2:leny-1,:)=( ...
        Array_water(:,3:leny,:)-2*Array_water(:,2:leny-1,:)+ ...
        Array_water(:,1:leny-2,:))/(dy*dy);
    A_h2o_i_ypp(:,2:leny-1,:)=( ...
        A_h2o_ion(:,3:leny,:)-2*A_h2o_ion(:,2:leny-1,:)+ ...
        A_h2o_ion(:,1:leny-2,:))/(dy*dy);
    %z的二阶偏导
    Array_water_zpp(:,:,2:lenz-1)=( ...
        Array_water(:,:,3:lenz)-2*Array_water(:,:,2:lenz-1)+ ...
        Array_water(:,:,1:lenz-2))/(dz*dz);
    A_h2o_i_zpp(:,:,2:lenz-1)=( ...
        A_h2o_ion(:,:,3:lenz)-2*A_h2o_ion(:,:,2:lenz-1)+ ...
        A_h2o_ion(:,:,1:lenz-2))/(dz*dz);
    %速度计算迭代
    A_w_v_z=A_w_v_z - gg*dt...
        + A_wa.*(0-A_w_v_z)*dt;
    %% 迭z计算水的扩散密度分布
    A_h2o_d=dt*(D_h2o.*( Array_water_xpp + ...     %扩散项
        Array_water_ypp+Array_water_zpp)*10^-10);
    A_h2o_v=dt*A_w_v_z.*Array_water_zp*dt*10^-3;   %垂直下降速度
    A_h2o_l=dt*k1*Array_water;                     %损失项
    A_h2o_r=Release(t0,dt,dx,dy,dz,Array_water,NvH2o*NA/dx/dy/dz/2*10^-15);%释放项
    Array_water_new=Array_water+A_h2o_r-A_h2o_l+A_h2o_v+A_h2o_d;%...
    %         + dt*(D_h2o.*( Array_water_xpp + ...
    %         Array_water_ypp+Array_water_zpp)*10^-10 ...
    %         +A_w_v_z.*Array_water_zp*dt*10^-3 ...
    %         -k1*Array_water ...
    %         + Release(t0,dt,dx,dy,dz,Array_water,NvH2o*NA/dx/dy/dz*10^-15)/dt);%
%     Array_water_new(:,:,lenz-1:lenz)=0; %边界条件设置
    for k=100:1:200
        if max(max(Array_water_new(:,:,k))-max(Array_water_new(:,:,k-1)))>0
            Array_water_new(:,:,k)=Array_water_new(:,:,k-1)/5.4366;
        end
    end
    Array_water_new(Array_water_new<0)=0;  %计算精度偏差 导致出现负数 矫正
    %% 计算电子密度变化
    %% h2o+离子密度
    A1=k1*Array_water.*A_o_ion*dt;%线性变化导致过衰减，设定单位时间长度内的最大衰减范围（指数衰减，这里可以偷懒）
    A2=k12.*A_h2o_ion.*A_e.*(300./A_e_T).^0.5*dt;
    if find(A1>=A_o_ion)
        A1(A1>A_o_ion)=0.99*A_o_ion(A1>A_o_ion);
    end
    if find(A2>=A_e)
        A2(A2>A_e)=0.98*A_e(A2>A_e);
    end
%     B1=dt.*cos(B_inc).*sin(B_inc).*cos(B_sit).*A_h2o_i_yp.*(D_h2o_ion./H_h2o_ion+2000);%*10^-5;   %随x方向的变量
    B1=dt.*B1_a.*A_h2o_i_yp.*(D_H ).*10^-5 ;   %随x方向的变量
%     B2=dt.*cos(B_inc).*sin(B_inc).*sin(B_sit).*A_h2o_i_xp.*(D_h2o_ion./H_h2o_ion+2000);%*10^-5 ;   %随y方向的变量
    B2=dt.*B2_a.*A_h2o_i_xp.*(D_H ).*10^-5;   %随y方向的变量
%         B3=dt*sin(B_inc)^2*(D_h2o_i_zp.*(A_h2o_i_zp +A_h2o_ion.*T_i_zp./T_i +A_h2o_ion./H_h2o_ion ) ...
%              +D_h2o_ion.*A_h2o_i_zp.*(T_i_zp./T_i + 0.01./H_h2o_ion) ... %单位补偿0.01./H
%              +D_h2o_ion.*A_h2o_ion.*(-1*(T_i_zp./T_i).^2 +T_i_zpp./T_i -H_h2o_i_zp./(H_h2o_ion.^2)*0.01) );
%     B3=dt.*sin(B_inc).^2.*(A_h2o_i_zp.*(D_h2o_i_zp + D_h2o_ion./T_i.*T_i_zp +D_h2o_ion./H_h2o_ion +2000));%*10^-5; %随z方向的变量
    B3=dt.*B3_a.*(A_h2o_i_zp.*(B3_a_1 +2000)).*10^-5; ...
%         +A_h2o_ion.*());%*10^-5; %随z方向的变量
%     B_D=dt.*(sin(B_inc).^2.*D_h2o_ion.*A_h2o_i_zpp ...   %扩散
%         +(cos(B_inc).*cos(B_sit)).^2.*D_h2o_ion.*A_h2o_i_ypp ...
%         +(cos(B_inc).*sin(B_sit)).^2.*D_h2o_ion.*A_h2o_i_xpp ).*10^-10;
    B_D=dt.*D_h2o_ion.*(B3_a.*A_h2o_i_zpp ...   %扩散
        +(B_i_cos.*B_d_cos).^2.*A_h2o_i_ypp ...
        +(B_i_cos.*B_d_sin).^2.*A_h2o_i_xpp ).*10^-10;

%     A_h2o_ion = A_h2o_ion ...
%         +A1 ... %第一步反应生成h2o++D_h2o_ion.*(A_h2oi_v_z)
%         -A2 ... %第二步损失项
%         +B_D ...
%         +B1+B2+B3  ;
    A_h2o_ion = A_h2o_ion +A1-A2  +B_D + B1+B2+B3 ;

    A_h2o_ion(A_h2o_ion<0)=0;%同上，计算精度误差导致负值
%     if find(A_h2o_ion<0)
%         disp(A_h2o_ion(A_h2o_ion<0));
%         disp('lizibudui');
%         break
%     end
    %% 电子密度计算
    A_e= A_e - A2 + dt*A_q - dt*A_o_ion.*B_o_los;%-B_D-B1-B2-B3;%第二步反应消耗电子
    A_o_ion=A_o_ion - A1 + dt*A_q - dt*A_o_ion.*B_o_los;  %O+氧原子离子的变化
    A_o_ion(A_o_ion<0)=0;
    A_e(A_e<0)=0;
    Array_water=Array_water_new;
    t0=t0+dt;%时间进度条
    if (t0)/tend-kkk>=0.01
        kkk=0.01+kkk;
        disp(t0);
        toc;
        disp("完成"+int2str(kkk*100)+'%');
        Cell_water{int8(kkk*100)}=gather(Array_water_new);
%         Cell_e{int8(kkk*100)}=gather(A_e);
%         Cell_o{int8(kkk*100)}=gather(A_o_ion);
        Cell_h2o_ion{int8(kkk*100)}=gather(A_h2o_ion);
        Cell_B1{int8(kkk*100)}=gather(B1);
        Cell_B2{int8(kkk*100)}=gather(B2);
        Cell_B3{int8(kkk*100)}=gather(B3);
        Cell_B_D{int8(kkk*100)}=gather(B_D);
%         filenameh2o=['.\mat\20241125\A_h2o_line',num2str(int8(tend/60)),'min',num2str(int8(kkk*100)),'.mat'];
%         filenamee=['.\mat\20241125\A_e_line',num2str(int8(tend/60)),'min',num2str(int8(kkk*100)),'.mat'];
%         filenameo=['.\mat\20241125\A_o_line',num2str(int8(tend/60)),'min',num2str(int8(kkk*100)),'.mat'];
%         filenameh2oi=['.\mat\20241125\A_h2oi_line',num2str(int8(tend/60)),'min',num2str(int8(kkk*100)),'.mat'];
%         A_h2o_i=gather(Array_water_new);
%         A_e_i=gather(A_e);
%         A_o_i=gather(A_o_ion);
%         A_h2oi_i=gather(A_h2o_ion);
%         save(filenameh2o,"A_h2o_i");
%         save(filenamee,"A_e_i");
%         save(filenameo,'A_o_i');
%         save(filenameh2oi,'A_h2oi_i');
%         clear A_h2o_i A_e_i A_o_i A_h2oi_i;
    elseif t0==dt
        disp('start');
    end
end


%% ATLAS-V火箭二级点火高度193km持续时间946-260s，火箭运动速度在4.7km/s以上持续加速
%火箭路径lon[-75=>-48],lat[28=>21]
%设置火箭释放尾气函数
function [ne]=Release(t,dt,dx,dy,dz,array,releaseN)
ne=0*array;
% xyz=int32(length(ne(1,1,:))/2);sa
if mod(t,0.5)<=dt  %每0.5秒释放一次
    lenx=length(array(:,1,1));
    leny=length(array(1,:,1));
    lenz=length(array(1,1,:));
    [x,y,z] = rocket_plume(lenx,leny,lenz,dx,dy,dz,t);
    x=round(x);y=round(y);z=round(z);
    if z>lenz||x>lenx||y>leny %超出范围就不释放了
        ne=0;
    else
        ne(x,y,z)=releaseN;
    end
else
    ne=0;
end
end

%%扩散系数：水：nO->氧原子浓度，nN2->氮气浓度，nO2->氧气浓度，T->背景温度
function [h2o]=Dh2o(nO,nN2,nO2,T)
h2o=(nO/(8.46*10^17 * T^0.5)+nN2/(2.04*10^17 * T^0.632) ...
    +nO2/(2.02*10^17 * T^0.632))^-1;
end
%扩散系数：水离子
function [h2o_i]=Dh2o_i(nO,nN2,nO2,T)
h2o_i=3.2*10^16.*T./(2.59*nO + 4.28*nO2 +4.39*nN2);
end
%扩散系数：O+
function [o_i]=Do_i(nO,nN2,nO2,T_i,T_n)
o_i=( (T_i+T_n).^0.5*nO./(3.1*10^17.*T_i) + nO2./(7.2*10^15.*T_i) ...
    +nN2./(8.5*10^15.*T_i) ).^-1;
end
%扩散系数：H2
function [H2]=Dh2(nO,nN2,nO2,T)
H2=(nO/(2.97*10^18 * T^0.5)+nN2/(2.8*10^17 * T^0.74) ...
    +nO2/(3.06*10^17 * T^0.732))^-1;
end
%太阳天顶角
function [szi]=solarze(lt,doy,lat)
pi=3.1415926536;
h=15*(lt-12)*pi/180;
delt=-23.44*cos(360*(doy+10)/365)/180*pi;
elev=asin(sin(delt)*sin(lat)+cos(delt)*cos(lat)*cos(h));
szi=acos((sin(delt)*cos(lat)-cos(delt)*cos(lat)*cos(h))/cos(elev))-pi/2;
end
%电子生成率
function [q]=epro(qm,hg,m,h,T,sza)
kk=1.380649*10^-23;
z=(h-hg).*m.*9.80000./kk./T;
q=qm.*(exp(1-z-sec(sza).*exp(-z)));
end

