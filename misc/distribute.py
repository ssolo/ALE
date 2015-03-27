from sys import *

CPUs=64
if len(argv)>2: CPUs=int(argv[2])
    
times={}
Cladess={}
Gss={}
Sss={}
fnames={}

for line in file(argv[1]).readlines():
    time,Clades,Gs,Ss,fname=line.strip().split()
    time=float(time)
    Clades=float(Clades)
    Gs=float(Gs)
    Ss=float(Ss)
    
    while times.has_key(time):
        time+=1e-6
    times[time]=time
    Gss[time]=Gs
    Sss[time]=Ss
    fnames[time]=fname
CPUs=min(CPUs,len(times.keys()))
bins={}
costs={}
for i in range(CPUs):
    bins[i]=[]
    costs[i]=0.

keeps=times.keys()



keeps=filter(lambda t: 0<Sss[t],keeps)
S=max(Sss.values())
near_universals=filter(lambda t: S*0.8<Sss[t],keeps)
i=0
nu_bound=0
for time in sorted(near_universals):
    i+=1
    if i/float(len(near_universals))>0.8:
        nu_bound=time
        break;
        
    
optimal_cost=max(times.keys())+1
keeps=filter(lambda t: optimal_cost>t,keeps)
optimal_cost=max([sum(keeps)/float(CPUs),nu_bound])
tmp=optimal_cost
print optimal_cost,max(keeps),-len(keeps)+len(times.keys())
while 1:
    tmp=optimal_cost
    old_keeps=keeps
    keeps=filter(lambda t: optimal_cost>t,keeps)
    optimal_cost=max([sum(keeps)/float(CPUs),nu_bound])
    print "optimal cost: ",optimal_cost,"max time: ",max(keeps)," discarding total of: ",-len(keeps)+len(times.keys()),"improvment: ",(tmp-optimal_cost)/tmp," discarded : ",int((1-(len(filter(lambda t: S*0.8<Sss[t],keeps))/float(len(near_universals))))*1000)/10.," % of ",len(near_universals)," near univ.s"

    if (tmp==optimal_cost or int((1-(len(filter(lambda t: S*0.8<Sss[t],keeps))/float(len(near_universals))))*1000)/10.>20.):
        keeps=old_keeps
        break
optimal_cost=max([sum(keeps)/float(CPUs),nu_bound])
keeps=filter(lambda t: optimal_cost>t,keeps)
print "stoped at ","optimal cost: ",optimal_cost,"max time: ",max(keeps)," discarding total of: ",-len(keeps)+len(times.keys()), " fams ",int((1-(len(filter(lambda t: S*0.8<Sss[t],keeps))/float(len(near_universals))))*1000)/10.," % of ",len(near_universals)," near univ.s"


bins={}
costs={}
for i in range(CPUs):
    bins[i]=[]
    costs[i]=0.
j=0
for time in sorted(keeps):
    bins[j%CPUs]+=[time]
    costs[j%CPUs]+=time
    j+=1

old_diff=-1
new_diff=max(costs.values())-min(costs.values())


while old_diff!=new_diff and max(costs.values())-min(costs.values())>1:
    bin_costs={}
    for i in range(CPUs):
        while bin_costs.has_key(costs[i]):
            costs[i]+=1e-6
        bin_costs[costs[i]]=i
    max_cost=sorted(bin_costs.keys())[-1]
    max_i=bin_costs[max_cost]

    for cost in sorted(bin_costs.keys()):
        other_i=bin_costs[cost]
        if other_i==max_i:
            pass;#stop=1
        
        for k in range(len(bins[max_i])):
            for j in range(len(bins[other_i])):
                max_time=bins[max_i][k]
                other_time=bins[other_i][j]
                if  abs( costs[max_i] - costs[other_i] ) > abs( (costs[max_i]-max_time+other_time) - (costs[other_i]+max_time-other_time) ) :
                    bins[max_i][k]=other_time
                    bins[other_i][j]=max_time
                    costs[max_i]-=max_time
                    costs[other_i]-=other_time
                    costs[max_i]+=other_time
                    costs[other_i]+=max_time
    old_diff=new_diff
    new_diff=max(costs.values())-min(costs.values())
    stdout.write("\r" "max: "+repr(max(costs.values()))[:5]+" min: "+repr(min(costs.values( )))[:5])
    stdout.flush()

print "achived max:",max(costs.values())," min:",min(costs.values())

fout=file(argv[1]+".conservative."+repr(CPUs),"w")
for i in range(CPUs):
    for time in bins[i]:
        fout.write(fnames[time]+"\t"+repr(i)+"\t"+repr(time)+"\n")
fout.close()


cut=1.
cut_keep=keeps
while len(cut_keep)>len(keeps)*0.8:
    cut-=0.001
    Sbins={}
    cut_keep=[]
    for time in keeps:
        if not Sbins.has_key(Sss[time]):
            Sbins[Sss[time]]=[]
        Sbins[Sss[time]]+=[time]
    for S in Sbins.keys():
        Ssum=sum(Sbins[S])
        cum_Ssum=0.
        for time in sorted(Sbins[S]):
            cum_Ssum+=time
            if cum_Ssum/Ssum<cut:
                cut_keep+=[time]
                
optimal_cost=max([sum(cut_keep)/float(CPUs),nu_bound])
print "cut with ",cut," discards ",-len(cut_keep)+len(keeps),"optimal cost: ",optimal_cost,"max time: ",max(cut_keep)," discarding total of: ",-len(cut_keep)+len(times.keys()), " fams ",int((1-(len(filter(lambda t: S*0.8<Sss[t],cut_keep))/float(len(near_universals))))*1000)/10.," % of ",len(near_universals)," near univ.s" 
keeps=cut_keep

bins={}
costs={}
for i in range(CPUs):
    bins[i]=[]
    costs[i]=0.
j=0
for time in sorted(keeps):
    bins[j%CPUs]+=[time]
    costs[j%CPUs]+=time
    j+=1

old_diff=-1
new_diff=max(costs.values())-min(costs.values())


while old_diff!=new_diff and max(costs.values())-min(costs.values())>1:
    bin_costs={}
    for i in range(CPUs):
        while bin_costs.has_key(costs[i]):
            costs[i]+=1e-6
        bin_costs[costs[i]]=i
    max_cost=sorted(bin_costs.keys())[-1]
    max_i=bin_costs[max_cost]

    for cost in sorted(bin_costs.keys()):
        other_i=bin_costs[cost]
        if other_i==max_i:
            pass;#stop=1
        
        for k in range(len(bins[max_i])):
            for j in range(len(bins[other_i])):
                max_time=bins[max_i][k]
                other_time=bins[other_i][j]
                if  abs( costs[max_i] - costs[other_i] ) > abs( (costs[max_i]-max_time+other_time) - (costs[other_i]+max_time-other_time) ) :
                    bins[max_i][k]=other_time
                    bins[other_i][j]=max_time
                    costs[max_i]-=max_time
                    costs[other_i]-=other_time
                    costs[max_i]+=other_time
                    costs[other_i]+=max_time
    old_diff=new_diff
    new_diff=max(costs.values())-min(costs.values())
    stdout.write("\r" "max: "+repr(max(costs.values()))[:5]+" min: "+repr(min(costs.values( )))[:5])
    stdout.flush()
print "achived max:",max(costs.values())," min:",min(costs.values())

fout=file(argv[1]+".aggressive."+repr(CPUs),"w")
for i in range(CPUs):
    for time in bins[i]:
        fout.write(fnames[time]+"\t"+repr(i)+"\t"+repr(time)+"\n")
fout.close()


#-----------------------------------


cut=1.
cut_keep=keeps
while len(cut_keep)>len(keeps)*0.25:
    cut-=0.001
    Sbins={}
    cut_keep=[]
    for time in keeps:
        if not Sbins.has_key(Sss[time]):
            Sbins[Sss[time]]=[]
        Sbins[Sss[time]]+=[time]
    for S in Sbins.keys():
        Ssum=sum(Sbins[S])
        cum_Ssum=0.
        for time in sorted(Sbins[S]):
            cum_Ssum+=time
            if cum_Ssum/Ssum<cut:
                cut_keep+=[time]
                
optimal_cost=max([sum(cut_keep)/float(CPUs),nu_bound])
print "cut with ",cut," discards ",-len(cut_keep)+len(keeps),"optimal cost: ",optimal_cost,"max time: ",max(cut_keep)," discarding total of: ",-len(cut_keep)+len(times.keys()), " fams ",int((1-(len(filter(lambda t: S*0.8<Sss[t],cut_keep))/float(len(near_universals))))*1000)/10.," % of ",len(near_universals)," near univ.s" 
keeps=cut_keep

bins={}
costs={}
for i in range(CPUs):
    bins[i]=[]
    costs[i]=0.
j=0
for time in sorted(keeps):
    bins[j%CPUs]+=[time]
    costs[j%CPUs]+=time
    j+=1

old_diff=-1
new_diff=max(costs.values())-min(costs.values())


while old_diff!=new_diff and max(costs.values())-min(costs.values())>1:
    bin_costs={}
    for i in range(CPUs):
        while bin_costs.has_key(costs[i]):
            costs[i]+=1e-6
        bin_costs[costs[i]]=i
    max_cost=sorted(bin_costs.keys())[-1]
    max_i=bin_costs[max_cost]

    for cost in sorted(bin_costs.keys()):
        other_i=bin_costs[cost]
        if other_i==max_i:
            pass;#stop=1
        
        for k in range(len(bins[max_i])):
            for j in range(len(bins[other_i])):
                max_time=bins[max_i][k]
                other_time=bins[other_i][j]
                if  abs( costs[max_i] - costs[other_i] ) > abs( (costs[max_i]-max_time+other_time) - (costs[other_i]+max_time-other_time) ) :
                    bins[max_i][k]=other_time
                    bins[other_i][j]=max_time
                    costs[max_i]-=max_time
                    costs[other_i]-=other_time
                    costs[max_i]+=other_time
                    costs[other_i]+=max_time
    old_diff=new_diff
    new_diff=max(costs.values())-min(costs.values())
    stdout.write("\r" "max: "+repr(max(costs.values()))[:5]+" min: "+repr(min(costs.values( )))[:5])
    stdout.flush()
print "achived max:",max(costs.values())," min:",min(costs.values())

fout=file(argv[1]+".skeleton."+repr(CPUs),"w")
for i in range(CPUs):
    for time in bins[i]:
        fout.write(fnames[time]+"\t"+repr(i)+"\t"+repr(time)+"\n")
fout.close()

