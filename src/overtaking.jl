type Car
    l::Float64
    acc::Float64
    v::Float64
    vmax::Float64
    pos::Float64
    over::Bool
end

type NS
    ρ::Float64
    p::Float64
    q::Float64
    L::Float64
    dt::Float64
    lsafe::Float64
end

function setup(;l=1.0, lsafe=7.0,acc=1.0, vmax=5.0, ρ=0.1, p=0.15, q=0.0, L=1000.0, dt=1.0)
    N = round(Int, ρ*L/(l+lsafe))
    NS(ρ,p,q,L,dt,lsafe), [Car(l,acc,0.0,vmax,pos,false) for pos in sort(sample(l:(l+lsafe):L,N,replace=false))]
end

function forward!(ns::NS, car::Car, observ::Float64)
    ispass = false
    if car.pos <= observ && car.pos + ns.dt*car.v > observ
        ispass = true
    end
    car.pos += ns.dt*car.v
    if car.pos > ns.L
        car.pos -= ns.L
    end
    ispass
end
function isin!(ns::NS, car::Car, observ1::Float64,observ2::Float64)
    isin = false
    if observ1 < car.pos <= observ2
        isin = true
    end
    isin
end

function onestep!(ns::NS, cars::Vector{Car}, cars1::Vector{Car}, ns_cars_idx::Vector{Int},
    nsos_cars_idx::Vector{Int}, idx::Vector{Int})
    overcnt = 0
    passcnt = 0
    incnt=0
    N = length(cars)
    resize!(ns_cars_idx,0)
    resize!(nsos_cars_idx,0)
    for i=2:N-1
        if rand() < ns.q
            push!(nsos_cars_idx, i)
        else
            push!(ns_cars_idx, i)
        end
    end
    append!(ns_cars_idx, [1;N])

    sort!(nsos_cars_idx, rev=true)

    for i in ns_cars_idx
        ip1 = i<N ? i+1 : 1
        ip2 = i<N-1 ? i+2 : (i<N ? 1 : 2)
        if cars[i].v + ns.dt*cars[i].acc > cars[i].vmax
            cars[i].v = cars[i].vmax
        else
            cars[i].v += ns.dt*cars[i].acc
        end
        d = cars[ip1].pos - cars[i].pos - cars[ip1].l-ns.lsafe
        if d<0
            d += ns.L
        end
        if d <= cars[i].v
            cars[i].v = d
        end
        if rand()<ns.p && cars[i].v >= ns.dt*cars[i].acc
            cars[i].v -= ns.dt*cars[i].acc
        end
    end

    for i in nsos_cars_idx
        ip1 = i<N ? i+1 : 1
        ip2 = i<N-1 ? i+2 : (i<N ? 1 : 2)
        ip3 = i<N-2 ? i+3 : (i<N-1 ? 1 : 2)
        if cars[i].v + ns.dt*cars[i].acc > cars[i].vmax
            cars[i].v = cars[i].vmax
        else
            cars[i].v += ns.dt*cars[i].acc
        end
        d = cars[ip1].pos - cars[i].pos - cars[ip1].l-ns.lsafe
        if cars[ip2].over
          d1 = cars[ip3].pos + ns.dt*cars[ip3].v - cars[ip1].pos - ns.dt*cars[ip1].v-cars[ip3].l-ns.lsafe
        else
          d1 = cars[ip2].pos + ns.dt*cars[ip2].v - cars[ip1].pos - ns.dt*cars[ip1].v-cars[ip2].l-ns.lsafe
        end
        if d<0
            d += ns.L
        end
        if d1<0
            d1 += ns.L
        end
        if cars[ip1].over
            cars[i].v = min(d + ns.dt*cars[ip1].v-cars[ip2].l, cars[i].v)
            if rand()<ns.p && cars[i].v >= ns.dt*cars[i].acc
                cars[i].v -= ns.dt*cars[i].acc
            end
        elseif d + cars[i].l + cars[ip1].l + ns.dt*cars[ip1].v <= ns.dt*cars[i].v && d1 >= cars[i].l
            cars[i].v = d + cars[i].l + cars[ip1].l + ns.dt*cars[ip1].v
            cars[i].over = true
            overcnt += 1
        else
            cars[i].v=min(d+ns.dt*cars[ip1].v,cars[i].v)
            if rand()<ns.p && cars[i].v >= ns.dt*cars[i].acc
                cars[i].v -= ns.dt*cars[i].acc
            end
        end
    end

    vsum = 0.0
    for i in eachindex(cars)
        if isin!(ns, cars[i], ns.L/2, ns.L/2+0.1*ns.L)
            incnt += cars[i].l
        end
        if forward!(ns, cars[i], ns.L/2)
            passcnt += 1
            vsum += cars[i].v
        end
    end
    for i in eachindex(cars)
        cars1[i] = cars[i]
    end
    carpos = [cars[i].pos for i in eachindex(cars)]
    idx = sortperm(carpos)
    for (i, pos) in enumerate(idx)
        cars[i] = cars1[pos]
        cars[i].over = false
    end
    overcnt, passcnt,incnt, vsum
end

function spacetime(file, ns, cars; pre_T=1000,T=100)
    ns_cars_idx = Int[]
    nsos_cars_idx = Int[]
    sizehint!(ns_cars_idx, length(cars))
    sizehint!(nsos_cars_idx, length(cars))
    idx = Array(Int, length(cars))
    cars1 = copy(cars)
    Road = fill(-1, T, round(Int,ns.L))
    for i=1:pre_T
        onestep!(ns, cars, cars1, ns_cars_idx, nsos_cars_idx, idx)
    end
    for t=1:T
        for car in cars
            Road[t, round(Int,car.pos)] = car.v
        end
        onestep!(ns, cars, cars1, ns_cars_idx, nsos_cars_idx, idx)
    end
    f = open(file, "w")
    for i=1:size(Road,1)
        for j=1:size(Road,2)
            v = Road[i,j] < 0 ? "*" : string(Road[i,j])
            write(f, v)
        end
        write(f, '\n')
    end
    close(f)
end
function trajectory(ns, cars; pre_T=1000,T=100)
    ns_cars_idx = Int[]
    nsos_cars_idx = Int[]
    sizehint!(ns_cars_idx, length(cars))
    sizehint!(nsos_cars_idx, length(cars))
    idx = Array(Int, length(cars))
    cars1 = copy(cars)
    Road = fill(-1, T, round(Int,ns.L))
    for i=1:pre_T
        onestep!(ns, cars, cars1, ns_cars_idx, nsos_cars_idx, idx)
    end
    for t=1:T
        for car in cars
          carend=car.pos-car.l+1
          if carend > 0
            Road[t, (round(Int,carend)):(round(Int,car.pos))] = car.v
          else
            Road[t,1:(round(Int,car.pos))] = car.v
            Road[t,(round(Int,carend+ns.L)):(round(Int,ns.L))] = car.v
          end
        end
        onestep!(ns, cars, cars1, ns_cars_idx, nsos_cars_idx, idx)
    end
    Road
end

function draw(ns::NS,car::Car, t=1, T=100)
    (context(),rectangle(car.pos/ns.L, t/T, car.l/ns.L, 1/ns.L),fill(RGB(1-car.v/car.vmax,0,0)))
end

function trajectory!(ns, cars; T=100)
    allcars = Any[]
    ns_cars_idx = Int[]
    nsos_cars_idx = Int[]
    sizehint!(ns_cars_idx, length(cars))
    sizehint!(nsos_cars_idx, length(cars))
    idx = Array(Int, length(cars))
    cars1 = copy(cars)
      for t=1:T
          for car in cars
              push!(allcars, draw(ns,car, t, T))
          end
          onestep!(ns, cars, cars1, ns_cars_idx, nsos_cars_idx, idx)
      end
      compose(context(), allcars...)
  end

function flowrate(ns, cars; pre_T=1000, T=2000)
    passcntsum = 0
    overcntsum = 0
    ns_cars_idx = Int[]
    nsos_cars_idx = Int[]
    sizehint!(ns_cars_idx, length(cars))
    sizehint!(nsos_cars_idx, length(cars))
    idx = Array(Int, length(cars))
    cars1 = copy(cars)
    for i=1:pre_T
        onestep!(ns, cars, cars1, ns_cars_idx, nsos_cars_idx, idx)
    end
    for i=1:T
        #a, b = onestep!(ns, cars, cars1, ns_cars_idx, nsos_cars_idx, idx)
        overcnt, passcnt,incnt, vsum=onestep!(ns, cars, cars1, ns_cars_idx, nsos_cars_idx, idx)
        overcntsum += overcnt
        passcntsum += passcnt
    end
    overcntsum/T/length(cars), passcntsum/T
end

function flux!(ns, cars; pre_T=1000, T=100*60)
    density = Float64[]
    flux = Float64[]
    velocity = Float64[]
    passcnt = 0
    overcnt = 0
    ns_cars_idx = Int[]
    nsos_cars_idx = Int[]
    sizehint!(ns_cars_idx, length(cars))
    sizehint!(nsos_cars_idx, length(cars))
    idx = Array(Int, length(cars))
    cars1 = copy(cars)
    for i=1:pre_T
        onestep!(ns, cars, cars1, ns_cars_idx, nsos_cars_idx, idx)
    end
    fluxsum = 0.0
    vsumsum = 0.0
    for i=1:T
        #a, b, c = onestep1!(ns, cars, cars1, ns_cars_idx, nsos_cars_idx, idx)
        overcnt, passcnt,incnt, vsum=onestep!(ns, cars, cars1, ns_cars_idx, nsos_cars_idx, idx)
        fluxsum += passcnt
        vsumsum += vsum
        if i%60 == 0
            push!(density,fluxsum*fluxsum/vsumsum)
            push!(flux,fluxsum)
            push!(velocity,vsumsum/fluxsum)
            fluxsum = 0.0
            vsumsum = 0.0
        end
    end
    density, flux, velocity
end

function timeheadway(ns, cars;pre_T=1000, T=2000)
    free=Float64[]
    sync=Float64[]
    jams=Float64[]
    passcnt = 0
    overcnt = 0
    N = length(cars)
    ns_cars_idx = Int[]
    nsos_cars_idx = Int[]
    sizehint!(ns_cars_idx, length(cars))
    sizehint!(nsos_cars_idx, length(cars))
    idx = Array(Int, length(cars))
    cars1 = copy(cars)
    for i=1:pre_T
        onestep!(ns, cars, cars1, ns_cars_idx, nsos_cars_idx, idx)
    end
    for t=1:T
        overcnt, passcnt,incnt, vsum=onestep!(ns, cars, cars1, ns_cars_idx, nsos_cars_idx, idx)
        density=incnt/(0.1*ns.L)
        for i in eachindex(cars)
            if isin!(ns, cars[i], ns.L/2,ns.L/2+50)
                ip1 = i<N ? i+1 : 1
                d = cars[ip1].pos - cars[i].pos
                if d<0
                    d += ns.L
                end
                if cars[i].v > 0
                    τ=d/cars[i].v
                    if 0<τ<10
                      density<=0.12 ? push!(free,τ) : (density >0.35 ? push!(jams,τ) : push!(sync,τ))
                      #ns.ρ<=0.1 ? push!(free,τ) : (ns.ρ >0.5 ? push!(jams,τ) : push!(sync,τ))
                    end

                end
            end
        end
    end
    free,sync,jams
end

#ns, cars = setup(p=0.25,q=0.5,ρ=0.6,vmax=5,l=1, L=1000)
#timeheadway("free1.csv","sync.csv1","jams.csv", ns, cars;pre_T=1000,T=10000)
#trajectory(ns, cars; T=100)
#ns, cars = setup(ρ=0.2, q=0.3, p=0.25, vmax=25, l=5)
#density,flux,velocity=flux!(ns, cars; pre_T=1000, T=100*60)
#using PyPlot
#plot(density,flux,"o")
