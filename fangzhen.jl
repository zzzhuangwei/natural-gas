include("BWRS.jl")
include("pipe.jl")
include("zs.jl")
include("yl.jl")
#*****************************************************条件*****************************************
x = [0.975 0 0.002 0 0.002 0 0 0 0 0 0 0 0 0 0 0.016 0.005 0]#气体组分
Qz = [4550*(10^4)/(60*60*24)]#体积流量m³/s
#管段：
Pz = [5e6]#压力Pa
Tz = [293.0]#温度K
Ke = 0.0177e-3#粗糙度m
D = 1.016#外径m
d = 0.9812#内径m
Th = 298#环境温度K
GDZC = [0 160e3 0 160e3 0 160e3 0 0 160e3 0 0]#管段长度m
L = 2000#步长m
crxs = 1.15#W/(㎡·K)
#压气站:

nmin = [6000 0 5160 0 5160 0 0 5160 0 0 5160]#转速可行区间r/min
nmax = [10500 0 9030 0 9030 0 0 9030 0 0 9030]
ε=0.001#精度
k=[4 0 2 0 2 0 0 2 0 0 2]#每个压气站的压缩机台数
n0=[8000 0 9030 0 9030 0 0 9030 0 0 9030]#额定转速
a1=[-1423.7 0 -243.41 0 -243.41 0 0 -243.41 0 0 -243.41]
b1=[4713 0 684.18 0 684.18 0 0 684.18 0 0 684.18]
c1=[6635.9 0 8677.2 0 8677.2 0 0 8677.2 0 0 8677.2]
a2=[-13.91 0 -2.38 0 -2.38 0 0 -2.38 0 0 -2.38]
b2=[68.03 0 18.96 0 18.96 0 0 18.96 0 0 18.96]
c2=[2.08 0 47.72 0 47.72 0 0 47.72 0 0 47.72]#多项式系数
n = []#转速
W = []#压缩机电机功率


#分输站进气口
Qjc = [0 0 0 0 0 0 126.51*(10^4)/(60*60*24) 0 0 126.51*(10^4)/(60*60*24) 0]

#反算
Qf = [4550*(10^4)/(60*60*24)]
Pf = [5e6]
Tf = [293.0]
GDZCf = [-80e3 0 -80e3 0 -100e3]
Qjcf = [0 506.04*(10^4)/(60*60*24) 0 -253.02*(10^4)/(60*60*24) 0 ]
Lf = -2000
#*****************************************仿真*********************
#反算last(Pf) last(Tf) last(Qf) GDZCf[j]
for j in 1:5
        if GDZCf[j] != 0
            (Pfs,Tfs) = pipe(last(Pf),last(Tf),last(Qf),Ke,D,d,Th,GDZCf[j],Lf,crxs,x)
            push!(Pf,Pfs)
            push!(Tf,Tfs)

        else
            Qfs = last(Qf)-Qjcf[j]
            push!(Qf,Qfs)
        end
end
Pd = [10e6 0 10e6 0 10e6 0 0 10e6 0 0 last(Pf)]#出口压力Pa
#正算
#压气站
for i in 1:11
        if k[i] != 0
                    (nc,Wc,Tc) = zs(last(Pz),Pd[i],last(Tz),last(Qz),n0[i],ε,nmin[i],nmax[i],x,a1[i],b1[i],c1[i],a2[i],b2[i],c2[i],k[i])
                            if nc > nmax[i]
                                (Pc1,Tc1,Wc1)=yl(last(Pz),last(Tz),nmax[i],last(Qz),k[i],x,a1[i],a2[i],b1[i],b2[i],c1[i],c2[i],n0[i])
                                    if Tc1>333
                                        Tc1=333
                                        push!(Pz,Pc1)
                                        push!(Tz,Tc1)
                                        push!(W,Wc1)
                                        push!(n,nmax[i])
                                    else
                                        push!(Pz,Pc1)
                                        push!(Tz,Tc1)
                                        push!(W,Wc1)
                                        push!(n,nmax[i])
                                    end
                            else
                                    if Tc>333
                                        Tc=333
                                        push!(Pz,Pd[i])
                                        push!(Tz,Tc)
                                        push!(n,nc)
                                        push!(W,Wc)
                                    else
                                        push!(Pz,Pd[i])
                                        push!(Tz,Tc)
                                        push!(n,nc)
                                        push!(W,Wc)
                                    end
                            end


        #管段
        elseif GDZC[i] != 0
            (Pc,Tc)=pipe(last(Pz),last(Tz),last(Qz),Ke,D,d,Th,GDZC[i],L,crxs,x)
            push!(Pz,Pc)
            push!(Tz,Tc)
        #分输站or进气口
        elseif Qjc[i] != 0
            Qc = last(Qz)-Qjc[i]
            push!(Qz,Qc)
        end


end

pop!(Pf)
pop!(Qf)
pop!(Tf)
Pfs = reverse(Pf)
Tfs = reverse(Tf)
Qfs = reverse(Qf)
append!(Pz,Pfs)
append!(Tz,Tfs)
append!(Qz,Qfs)
