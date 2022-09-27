#RISM solver using Picard iteration scheme for hardsphere dimers

using FFTW
using Plots
using LinearAlgebra

function main()


    #define variables
    #define grid size with N
    #define grid spacing with delr
    N = 2048
    delr = 0.04
    delk = pi/delr/N
    alpha = 0.1
    converge = 1.0*10^(-7) 
    den = 0.0191 
    count = 0
    
    #define dimer site sizes and bondlength
    sigma_A = 0.0 
    sigma_B = 3.47
    sigma_AB = 1.17
    bond_L = 2.3
    #(from above this is signifying a tangential bond)
    # r and k values
    rvalues = ((1:N).-0.5)*delr
    kvalues = ((1:N).-0.5)*delk
    #println(rvalues)
    #println(kvalues)
    #hrclosure vecotrs
    AA_close_index = round(Int,sigma_A/delr)
    BB_close_index = round(Int,sigma_B/delr)
    AB_close_index = round(Int,sigma_AB/delr)
    BA_close_index = AB_close_index
    #println(AB_close_index)
    
    omegakAA = ones(N)
    omegakBB = omegakAA
    omegakAB = sin.(kvalues*bond_L)./(kvalues*bond_L)
    omegakBA = sin.(kvalues*bond_L)./(kvalues*bond_L)
    #println(omegak)
    #a = plot(kvalues,omegak,seriestype =:scatter, title = "testing omega_k plot",label = "f = sin(kL)/kL",legend =:top)
    #display(a)
    #readline allos for the plot to be seen
    #readline()
    #println(omegaAA)
    #let's creat the vector of omega matrix elements from the above vectors
    omegak = zeros(2,2,N)
    for i in 1:N
        omegak[1,1,i] = omegakAA[i]
        omegak[1,2,i] = omegakAB[i]
        omegak[2,1,i] = omegakBA[i]
        omegak[2,2,i] = omegakBB[i]
    end
    #println(omegak[:,:,1])
    gammar = zeros(2,2,N)
    gammar_new = zeros(2,2,N)
    gammak = zeros(2,2,N)
    cr = zeros(2,2,N)    
    ck = zeros(2,2,N)    
    hr = zeros(2,2,N)
    hk = zeros(2,2,N)
    id = Matrix(1.0I,2,2)
    rho = id.*den
    inverse = zeros(2,2,N)
    AFD = 1.0

    #start of solving loop 
    #for z in 1:1000
        #find c(r) using the identity c(r) = h(r) - gamma(r) and apply the closures for c(r) and h(r)
    while AFD > converge 
        hr[1,1,1:AA_close_index].= -1.0
        hr[1,2,1:AB_close_index].= -1.0
        hr[2,1,1:BA_close_index].= -1.0
        hr[2,2,1:BB_close_index].= -1.0
        cr = hr - gammar
        cr[1,1,AA_close_index + 1: N].= 0.0
        cr[1,2,AB_close_index + 1: N].= 0.0
        cr[2,1,BA_close_index + 1: N].= 0.0
        cr[2,2,BB_close_index + 1: N].= 0.0
        #println(cr[1,1,:])
        #multiply by rvalues and divide by k values for each matrix element
        cr[1,1,:].*=rvalues
        cr[1,2,:].*=rvalues
        cr[2,1,:].*=rvalues
        cr[2,2,:].*=rvalues
        #println(cr[1,2,:])
        #perform the FFT on these data
        ck[1,1,:]=(FFTW.r2r(cr[1,1,:],FFTW.RODFT11)./kvalues)*2*pi*delr
        ck[1,2,:]=(FFTW.r2r(cr[1,2,:],FFTW.RODFT11)./kvalues)*2*pi*delr
        ck[2,1,:]=(FFTW.r2r(cr[2,1,:],FFTW.RODFT11)./kvalues)*2*pi*delr
        ck[2,2,:]=(FFTW.r2r(cr[2,2,:],FFTW.RODFT11)./kvalues)*2*pi*delr
        #println(ck[1,1,1])
        #println(ck[1,2,1])
        #println(ck[2,1,1])
        #println(ck[2,2,1])
        # this will look weird but i am just trying to copy my mathematica code...i will sort out the details later
        for d in 1:N
            gammak[:,:,d]= (omegak[:,:,d]*ck[:,:,d]*inv(id-rho*omegak[:,:,d]*ck[:,:,d])*omegak[:,:,d]-ck[:,:,d])*kvalues[d]
        end
        #println(gammak[:,:,1])
        gammar_new[1,1,:]=FFTW.r2r(gammak[1,1,:],FFTW.RODFT11)./(rvalues*4*N*pi*delr)
        gammar_new[1,2,:]=FFTW.r2r(gammak[1,2,:],FFTW.RODFT11)./(rvalues*4*N*pi*delr)
        gammar_new[2,1,:]=FFTW.r2r(gammak[2,1,:],FFTW.RODFT11)./(rvalues*4*N*pi*delr)
        gammar_new[2,2,:]=FFTW.r2r(gammak[2,2,:],FFTW.RODFT11)./(rvalues*4*N*pi*delr)
        #println(gammar_new[:,:,1])
        AFD = sqrt(1/N*((sum(((gammar[1,1,:].-gammar_new[1,1,:])./(gammar[1,1,:].+gammar_new[1,1,:])).^2)+sum(((gammar[1,2,:].-gammar_new[1,2,:])./(gammar[1,2,:].+gammar_new[1,2,:])).^2)+sum(((gammar[2,1,:].-gammar_new[2,1,:])./(gammar[2,1,:].+gammar_new[2,1,:])).^2)+sum(((gammar[2,2,:].-gammar_new[2,2,:])./(gammar[2,2,:].+gammar_new[2,2,:])).^2))))
        count +=1
        #thing = sum(((gammar[1,1,:]-gammar_new[1,1,:])/(gammar[1,1,:]+gammar_new[1,1,:]))^2)
        #:println(thing)
        #println(gammar[:,:,1])
        #println(gammar_new[:,:,1])
        gammar = gammar.*(1.0-alpha) +gammar_new.*(alpha)
        #println(gammar[:,:,1])
        #println(gammar[:,:,1])
        #println(AFD)
    end  
    #println(AA_close_index)
    #println(sigma_A)
    #println(sigma_A/delr)
    hr = gammar.+cr
    hr[1,1,1:AA_close_index].= -1.0
    hr[1,2,1:AB_close_index].= -1.0
    hr[2,1,1:BA_close_index].= -1.0
    hr[2,2,1:BB_close_index].= -1.0
    gr = hr .+ 1.0
    println(AFD)
    println(count)
    a = plot(rvalues,gr[1,1,:],seriestype=:scatter,xlims=(0,12),ylims=(0.0,3.5))
    b = plot(rvalues,gr[1,2,:],seriestype=:scatter,xlims=(0,12),ylims=(0.0,3.5))
    c = plot(rvalues,gr[2,2,:],seriestype=:scatter,xlims=(0,12),ylims=(0.0,3.5))
    d = plot(a,b,c,lagout = (3,1))
    display(d)
    readline()
end

main()
