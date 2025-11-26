using ReachabilityAnalysis, Plots, LazySets
using LaTeXStrings
plotlyjs(size = (1000, 1000))


@taylorize function duffing!(dx,x, params, t)
    local u1, u2, alph, alpha_e,alpha_h, T_e, T_h=0.5, 0.5, 0.01, 0.04, 0.145, 15.0, 40.0
    v0 = u1
    v1 = u2
    x0=x[1]
    x1=x[2]
    a_00 = 1 - 2 * alph - alpha_e - alpha_h * v0
    a_11 = 1 - 2 * alph - alpha_e - alpha_h * v1
    dx[1] = a_00 * x0 + alph * x1 + alpha_h*T_h*u1 + alpha_e*T_e
    dx[2] = a_11 * x1 + alph * x0 + alpha_h*T_h*u2 + alpha_e*T_e
    return dx
end

function sys(x,y)
    local u1, u2, alph, alpha_e,alpha_h, T_e, T_h=0.5, 0.5, 0.01, 0.04, 0.145, 15.0, 40.0
    v0 = u1
    v1 = u2
    x0=x
    x1=y
    a_00 = 1 - 2 * alph - alpha_e - alpha_h * v0
    a_11 = 1 - 2 * alph - alpha_e - alpha_h * v1
    temp1 = a_00 * x0 + alph * x1 + alpha_h*T_h*v0 + alpha_e*T_e
    temp2 = a_11 * x1 + alph * x0 + alpha_h*T_h*v1 + alpha_e*T_e
    return [temp1,temp2]
end

function rect_plot(H::Hyperrectangle,zaxi, col::String)
    c=LazySets.center(H)
    r=LazySets.radius_hyperrectangle(H)
    x = [c[1]-r[1], c[1]-r[1], c[1]+r[1], c[1]+r[1],
     c[1]-r[1], c[1]-r[1], c[1]+r[1], c[1]+r[1]]
    y = [c[2]-r[2], c[2]+r[2], c[2]-r[2], c[2]+r[2],
    c[2]-r[2], c[2]+r[2], c[2]-r[2], c[2]+r[2]]
    z = [0, 0, 0, 0, 
    zaxi, zaxi, zaxi, zaxi]

    connections = [(1,2,3),(2,3,4),(5,6,7),(6,7,8), 
                (2,4,8),(2,8,6),(1,2,6),(1,6,5),
                (3,4,7),(4,7,8),(1,3,7),(1,5,7)];


    mesh3d!(x,z,y; connections, linewidth=0,fill = col, label="", margin = 2 * Plots.mm, fillalpha=3)
end

function rect_plotz(H::Hyperrectangle,i, col::String, alp)
    
    c=LazySets.center(H)
    r=LazySets.radius_hyperrectangle(H)
    x = [c[1]-r[1], c[1]-r[1], c[1]+r[1], c[1]+r[1],
     c[1]-r[1], c[1]-r[1], c[1]+r[1], c[1]+r[1]]
    y = [c[2]-r[2], c[2]+r[2], c[2]-r[2], c[2]+r[2],
    c[2]-r[2], c[2]+r[2], c[2]-r[2], c[2]+r[2]]
    z = [i, i, i, i, 
    i, i, i, i]

    connections = [(1,2,3),(2,3,4),(5,6,7),(6,7,8), 
                (2,4,8),(2,8,6),(1,2,6),(1,6,5),
                (3,4,7),(4,7,8),(1,3,7),(1,5,7)];

    mesh3d!(x,z,y; connections, linewidth=0,fill = col, margin = 2 * Plots.mm, fillalpha=alp, label="")
end

function rect_plotem(H::Hyperrectangle,i, col::String, alp)
    local long=0.3
    
    c=LazySets.center(H)
    r=LazySets.radius_hyperrectangle(H)
    x = [c[1]-r[1], c[1]-r[1], c[1]+r[1], c[1]+r[1],
     c[1]-r[1], c[1]-r[1], c[1]+r[1], c[1]+r[1]]
    y = [c[2]-r[2], c[2]+r[2], c[2]-r[2], c[2]+r[2],
    c[2]-r[2], c[2]+r[2], c[2]-r[2], c[2]+r[2]]
    z = [i-long, i-long, i-long, i-long, 
    i+long, i+long, i+long, i+long]

    connections = [(1,2,3),(2,3,4),(5,6,7),(6,7,8), 
                (2,4,8),(2,8,6),(1,2,6),(1,6,5),
                (3,4,7),(4,7,8),(1,3,7),(1,5,7)];

    mesh3d!(x,z,y; connections, linewidth=0,fill = col, fillalpha=alp,label="")
end


const inter=16

z=collect(0:1:inter-1)

x_s=Array{Float64}(undef, inter) 
y_s=Array{Float64}(undef, inter) 
x_s[1]=20
y_s[1]=20
x=[x_s[1],y_s[1]]

delta=[0.5,0.5]

F_c=Float64[24, 24]
F_r=Float64[26, 26]

XF = Hyperrectangle(low=F_c, high=F_r)


#system trace
for i=2:inter
    #temp=f(x)
    local temp
    temp=sys(x_s[i-1],y_s[i-1])
    x_s[i]=temp[1]
    y_s[i]=temp[2]
    if temp ∈ XF
        print("Enter Fault region at step: ", i, "\n")
    end
    i=i+1  
end

#field of imprecise observation
Y_s=Vector{Hyperrectangle}(undef, inter) 

for i=1:inter
    center=Float64[x_s[i],y_s[i]]
    temp=Hyperrectangle(center, delta)
    
    Y_s[i]=temp
    i=i+1
end


#sample the imprecise observation from above field
x_i=Array{Float64}(undef, inter) 
y_i=Array{Float64}(undef, inter) 
for i=1:inter
    temp=sample(Y_s[i], 1)[1]
    x_i[i]=temp[1]
    y_i[i]=temp[2]
end

#field of imprecise observation
Y_i=Vector{Hyperrectangle}(undef, inter) 

for i=1:inter
    center=Float64[x_i[i],y_i[i]] 
    temp=Hyperrectangle(center, delta) 
    Y_i[i]=temp
    i=i+1
end

#Overall record the observation of Diagnoser, its index is time step
Re=[]

#step of fault alarm
index_f=0

# record the compute time
#
time_reach=0
# time_var = @elapsed begin



#Compute M(0)
m=difference(Y_i[1],XF) #UnionSetArray
#current observation region
temp=array(m) #Array
#push the initial region in overall Record

push!(Re,temp) # [ [temp1], [....] ]


#Compute the reachable record M1(0) 
M1=[] # vector of arrays
push!(M1,Re[1]) # 

M_c=temp
# print(length(M_c))
time_var = @elapsed begin
for i=2:inter
    # print(length(M1))
    M1_n=Hyperrectangle[]
    M_n=Hyperrectangle[]
    global M_c, Re,M1, index_f, time_reach
    # print(M_c)
    print(length(M_c),"\n")
    if length(M_c)==0
        break
    end
    for j =1:length(M_c) 
        #extract a pice of current region
        X_c=M_c[j] # set it as the initial region
        if X_c==∅
            continue
        end

        time_reach_t = @elapsed begin

        prob = @ivp(X' = duffing!(X), X(0) ∈ X_c, dim=2)
        sol = solve(prob, T=0.05, alg=TMJets21a()) # Reachsolution
        R = sol[1] # TaylorModel Reachset
        temp=overapproximate(R,Hyperrectangle) # Reachset 
        temp4=set(temp) # a Hyperrectangle
        #record element of M1_n
        if temp4==∅
            print("t4 is empty. \n")
            continue
        end
        push!(M1_n,temp4)

        end

        print("time-reach_t:", time_reach_t,"\n")
        time_reach=time_reach+time_reach_t

        #with fault
        temp2=intersection(temp4,Y_i[i]) # a Hyperrectangle which can be reached by normal behaviors
        
        if temp2==∅(2)
            print("t2 is empty. \n")
            continue
        end

        #eliminate fault
        temp3= difference(temp2,XF) # UnionSetArray
        if length(temp3)==0
            print("t3=empty \n")
            continue
        end
        # translate UnionSetArray into array to record and opearte
        temp5=array(temp3) # Hyperrectangle Array

        # record the next observation of the current region j in the cash M_next
        append!(M_n,temp5) #Hyperrectangle[]
        if length(M_n)>1
            A=[]
        end
    end

    # Record the next observation
    if length(M_n)==0
        print("\n M_n=empty in ", i, "th inter.\n")
        index_f=i-1 # record the fault step
        break
    else
        push!(Re,M_n)
        push!(M1,M1_n)
        M_c=M_n
    end
    
end

end

time_aver=time_var/index_f
print("time_aver: ", time_aver, "\n")

print("time_reach: ", time_reach, "\n")

time_reach_aver=time_reach/index_f
print("time_reach_aver: ", time_reach_aver, "\n")

########################################################
#Plots

# parent = Scene()
# cam3d!(parent)

# # One can set the camera lookat and eyeposition, by getting the camera controls and using `update_cam!`
# camc = cameracontrols(parent)
# update_cam!(parent, camc, Vec3f(0, 8, 0), Vec3f(4.0, 0, 0))

# plot(x_i,z,y_i, seriestype=:scatter,xlabel="x",ylabel="z",zlabel="y",markersize = 2, linewidth=3 ,mc=:red, legend=false)
# default(legendfontsize = 15, guidefont =18, tickfont = 10)
pl=plot(x_s,z,y_s; grid=true, seriestype=:path3d, linestyle=:solid,linewidth=2, linecolor =:green, markershape = :circle,
 markersize = 1.5,  mc=:green, markerstrokewidth=0,label="System's state run",legendfontsize = 15, guidefont =18, tickfont = 10, xlabel=L"x_1",ylabel=L"t",zlabel=L"x_2")    

plot!(x_i,z,y_i,xlims=(19,28),ylims=(0,inter),zlims=(19,28),xticks=19:1:28,zticks=19:1:28,yticks=0:1:inter, font=(15, "Courier"); seriestype=:scatter, 
    markershape = :circle, markersize = 1.5,  mc=:red, markerstrokewidth=0,label="Imprecise observation",legend=:top)

    # scatter!(x_i,z,y_i; yticks=0:1:inter,xlabel="x1",ylabel="t",zlabel="x2", legend=:topleft,"imprecise observation" , lab = "system state"
    # markershape = :circle, markersize = 1.5,  mc=:red, legend=false, draw_arrow = true,markerstrokewidth=0)

# plot XF
for i=1:inter
    rect_plotz(XF,i-1, "gray50",0.3)
end

for i=1:length(Re)
    for j in Re[i]
            rect_plotz(j,i-1, "blue4",0.6)
    end
    
    if i==length(Re) && i<inter
        local c, d =[x_i[i+1],y_i[i+1]], [0.3,0.3]
        em=Hyperrectangle(c,d)
        rect_plotem(em,i,"yellow",0.6)
        # l1=[x_i[i+1],x_i[i+1],x_i[i+1]+100]
        # l2=[i,i,i]
        # l3=[y_i[i+1],20,20]
        l1=[x_i[i+1],x_i[i+1]]
        l2=[i,i]
        l3=[y_i[i+1],15]
        plot!(l1,l2, l3; seriestype=:path, linestyle=:solid,linewidth=2, linecolor =:black,label="")
    end 
end

# show(pl)
display(pl)
gui()
# savefig("myplot.png") 

readline()
