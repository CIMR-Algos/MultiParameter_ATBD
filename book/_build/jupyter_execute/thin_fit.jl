using Pkg
Pkg.activate("../mpr_atbd_env_jl")
using CSV
using DataFrames
using ZipFile
using Dates
using LsqFit
using ForwardDiff
using StaticArrays
using YAML
#using PythonCall
using PythonPlot
#mpl=pyimport("matplotlib")
#mpl.use("svg")
#plt=pyimport("matplotlib.pyplot")
plt=PythonPlot.pyplot

I(p)=f(x)=p[2]-(p[2]-p[1])*exp(-x/p[3])
Q(p)=f(x)=p[2]-(p[2]-p[1])*exp(-(x/p[3])^p[4])


#reading fit parameters from file

fI(x)=I(pI)(x)
fQ(x)=Q(pQ)(x)

params=YAML.load_file("../../SeaIceThickness_ATBD/book/fit_params.yml")
const tbh_params,tbv_params,pI,pQ=getindex.(Ref(params),["ph","pv","pI","pQ"])


ff(p)=return f(x)=p[2]-(p[2]-p[1])*exp(-x/p[3])
sit_h=ff(tbh_params)
sit_v=ff(tbv_params)

Fw_TB(x,p)=SA[ff(p[1])(x), ff(p[2])(x)]
Fw_IQ(x)=[fI(x[1]), fQ(x[1])]

Fw_TB(x)=Fw_TB(x[1],(tbh_params,tbv_params))

function lm_retrieval(Ta,Sₑ,Sₐ,xₐ,F)
    #Levenberg Marquardt method after Rodgers (2000)
    #target: find x so that F(x)=Ta, given
    #Ta: measurement vector
    #Sₑ: error covariance of measurement
    #Sₐ: error covariance of physical state 
    #xₐ: expected physical state (also used as start, i.e. first guess)
    #F: the forward model translating measument space into state space
    Sₐ⁻¹=inv(Sₐ)
    Sₑ⁻¹=inv(Sₑ)
    #function to minimize with changing input x
    J(y,x,Sₑ⁻¹,Sₐ⁻¹,xₐ,F)=(y.-F(x))'*(Sₑ⁻¹*(y.-F(x)))+(xₐ.-x)'*(Sₐ⁻¹*(xₐ.-x)) 
    xᵢ=copy(xₐ)
    Jᵢ=J(Ta,xᵢ,Sₑ⁻¹,Sₐ⁻¹,xₐ,F)
    γ=1e-5 #set to 0 for gauss newton
    for i=1:2000 
        Kᵢ=ForwardDiff.jacobian(F,xᵢ)
        Ŝ⁻¹=Sₐ⁻¹+Kᵢ'*Sₑ⁻¹*Kᵢ #eq 5.13
        xᵢ₊₁=xᵢ+((1+γ)*Sₐ⁻¹+Kᵢ'*Sₑ⁻¹*Kᵢ)\(Kᵢ'*Sₑ⁻¹*(Ta-F(xᵢ))-Sₐ⁻¹*(xᵢ-xₐ)) #eq 5.36
        Jᵢ₊₁=J(Ta,xᵢ₊₁,Sₑ⁻¹,Sₐ⁻¹,xₐ,F)
        d²=(xᵢ-xᵢ₊₁)'*Ŝ⁻¹*(xᵢ-xᵢ₊₁) #eq 5.29
        if Jᵢ₊₁<Jᵢ 
            γ/=2
        else
            γ*=10
            continue
        end
        xᵢ=xᵢ₊₁
        if d²<1e-10
            break
        end
        Jᵢ=Jᵢ₊₁
    end
    Kᵢ=ForwardDiff.jacobian(F,xᵢ)
    Ŝ=inv(Sₐ⁻¹+Kᵢ'*Sₑ⁻¹*Kᵢ) # eq 5.38
    
    return xᵢ,Ŝ
end


retrievallm(h,v)=first.(lm_retrieval(SA[h,v],SA[25 15;15 25.0],SMatrix{1,1,Float64,1}(20000.0),SA[100.],Fw_TB))



zfn="../data/RRDP_v3.zip"
if !isfile(zfn)
    download("https://figshare.com/ndownloader/files/31422043",zfn)
end

zf=ZipFile.Reader(zfn);
idxs=findall(x->(occursin("RRDP_v3/UB_THINICE/S",x)),getfield.(zf.files,:name))
D=DataFrame[]
for idx in idxs
    push!(D,CSV.read(read(zf.files[idx]),DataFrame;skipto=3,header=2,missingstring=["-999","noval"]))
#    fio,DataFrame;skipto=3,header=2)
end

allD=vcat(D...);
sort!(allD,["date"]);

d=allD.date[1]
datef=dateformat"yyyymmdd"
Dates.format(d,dateformat"yyyymmdd")

for r in eachrow(allD)
    cday=Dates.format(r.date,datef)
  #  println("$(r.latitude) $(r.longitude) 51/55 $cday $cday")
end
rrdp_extract=CSV.read("../data/rrdp_thinice_51-55_724.txt",DataFrame,header=["date","tbh53","tbv53","mean_inc","tbh53_std","tbv53_std"]);
fulltable=hcat(allD,rrdp_extract,makeunique=true);
replace!( x->ismissing(x) ? NaN : x,fulltable[!,"89.0H"]) .|> Float64
replace!( x->ismissing(x) ? NaN : x,fulltable[!,"89.0V"]) .|> Float64
replace!(x->ismissing(x) ? NaN : x,fulltable[!,"SIT"]) .|> Float64
replace!(x->ismissing(x) ? NaN : x,fulltable[!,"sigma_40"]) .|> Float64
fulltable[!,"newSIT"]=retrievallm.(fulltable[!,"tbh53"],fulltable[!,"tbv53"]) .|> first
fulltable[!,"newSIT_std"]=retrievallm.(fulltable[!,"tbh53"],fulltable[!,"tbv53"]) .|> last
iidx=isnan.(fulltable[!,"tbh53"])
#fulltable[iidx,"newSIT"].=NaN
#fulltable[iidx,"newSIT_std"].=NaN
fulltable=fulltable[.!iidx,:]

#iidx=fulltable[!,:SIT].!=-999
#fulltable=fulltable[iidx,:];

retrievallm.(fulltable[!,"tbh53"],fulltable[!,"tbv53"])

fulltable[!,"tbh53"]

print(names(fulltable))

include("../algorithm/src/OEM_with_SIT.jl")
using .OEM

fig,axs=plt.subplots(nrows=7,ncols=2,figsize=(12,10))
#scatter(fulltable[!,:SIT],fulltable[!,"10.7H"])
channels=["6.9V","6.9H","10.7V","10.7H","18.7V","18.7H","23.8V","23.8H","36.5V","36.5H","89.0V","89.0H"]
owtp=[157.9397756705356
  74.22054192317813
 163.4896743050617
  78.40477088765373
 175.18191921420814
  90.60124656474672
 185.91398347505353
  99.89685544430296
 206.667948239106
 128.36037504633393
 242.01900700540162
 175.158277073163]
for i=1:12
    nanidx=.!isnan.(fulltable[!,channels[i]])
    axs.flat[i-1].scatter(fulltable[nanidx,:newSIT],fulltable[nanidx,channels[i]],c=fulltable[nanidx,:strd],cmap="viridis")
    cf(x,pp)=ff([owtp[i],pp[1],pp[2]]).(x)
    hp=curve_fit(cf,fulltable[nanidx,:newSIT],fulltable[nanidx,channels[i]],[200,2.0]).param
    hp=[owtp[i],hp...]
    #println(hp[3])
    #axs[p].scatter(fulltable[!,:newSIT].* 100,fulltable[!,channels[i]])
    xx=collect(0:1:150)
    axs.flat[i-1].plot(xx,ff(hp).(xx),color="red")
    #plt.colorbar(axs.flat[i])
end
gcf()

ax=axs.flat[1]
ax.


fig,ax=plt.subplots(figsize=(8,8))
ax.plot(fulltable[!,:SIT],".")
display(MIME("image/png"),fig)

fulltable[!,"89.0V"]



fn="../data/THINICE_SIC1-AMSR-N_wsmos.text"
mdf=CSV.read(fn,DataFrame)
sit=mdf[!,6]
for i=1:43
    println("$i $(names(mdf)[i])")
end

@. ff(x,p)=p[2]-(p[2]-p[1])*exp(-x/p[3])
sit=mdf[!,6]

idx= .!(sit==-999) .& .!(sit.>0.5)



fig,axs=plt.subplots(figsize=(20,5),nrows=2,ncols=6)
axs=pyconvert(Array,axs)
owtp=[157.9397756705356
  74.22054192317813
 163.4896743050617
  78.40477088765373
 175.18191921420814
  90.60124656474672
 185.91398347505353
  99.89685544430296
 206.667948239106
 128.36037504633393
 242.01900700540162
 175.158277073163] #ow tiepoints
for p=1:12
    #println(p)
    index=11+p
#idexes=
    dv=mdf[!,index]
    name=names(mdf)[index]
    cf(x,pp)=ff(x,[owtp[p],pp[1],pp[2]])
    hp=curve_fit(cf,sit.*100,dv,[200,2.0]).param
    hp=[owtp[p],hp...]
    #println(hp[3])
    axs[p].scatter(sit.* 100,dv)
    xx=collect(0:1:30)
    axs[p].plot(xx,ff.(xx,Ref(hp)),color="red")
    axs[p].set_title(name)
    axs[p].set_ylim(70,300)
    axs[p].grid()
end
display(MIME("image/png"),fig)

fulltable.SIT


