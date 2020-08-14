using DifferentialEquations, Plots

#Reference: https://diego.assencio.com/?index=e5ac36fcb129ce95a61f8e8ce0572dbf


#Solving the double pendulum with a traditional ODE method
#==========================================================#
function doublependulum(du, u, params, t)
    l1 = params[1]
    l2 = params[2]
    m1 = params[3]
    m2 = params[4]
    g  = params[5]

    p1 = u[1]
    p2 = u[2]
    θ1 = u[3]
    θ2 = u[4]

    h1 = p1*p2*sin(θ1-θ2)/(l1*l2*(m1+m2*sin(θ1-θ2)^2))
    h2 = (m2*l2^2*p1^2 + (m1 + m2) * l1^2 * p2^2 - 2*m2*l1*l2*p1*p2*cos(θ1-θ2))/(2 * l1^2 * l2^2 *(m1 + m2*sin(θ1-θ2)^2)^2)

    d_p1 = -(m1 + m2)* g*l1*sin(θ1) - h1 + h2*sin(2*(θ1-θ2))
    d_θ1 = (l2*p1 - l1*p2*cos(θ1-θ2))/(l1^2 * l2*(m1 + m2*sin(θ1 - θ2)^2))
    d_p2 = -m2*g*l2*sin(θ2) + h1 - h2*sin(2*(θ1-θ2))
    d_θ2 = (-m2*l2*p1*cos(θ1 - θ2) + (m1 + m2)*l1*p2) / (m2*l1*l2^2*(m1 + m2*sin(θ1-θ2)^2))


    du .= [d_p1, d_p2, d_θ1, d_θ2]
    return nothing
end

l1 = 1. #length of pendulum1
l2 = 2. #length of pendulum2
m1 = 1. #mass of pendulum1
m2 = 1. #mass of pendulum2
g = 9.81 #gravity

params = [l1, l2, m1, m2, g]
times = (0., 25.)
u0 = [1.,1.,1.,1.]
prob = ODEProblem(doublependulum, u0, times, params)
sol1 = solve(prob, AutoVern7(Rodas5()), dt = .005)

#plot solution
plot(sol1, vars=1, xlim=(0,20), label="Momentum1")
plot!(sol1, vars=2, xlim=(0,20), label="Momentum2")
plot!(sol1, vars=3, xlim=(0,20), label="theta1")
plot!(sol1, vars=4, xlim=(0,20), label="theta2")
#==========================================================#


#Now with HamiltonianProblem()
#==========================================================#
function H(p, θ, params)
    l1 = params[1]
    l2 = params[2]
    m1 = params[3]
    m2 = params[4]
    g  = params[5]

    return  (m2*l2^2*p[1]^2 + (m1+m2)*l1^2*p[2]^2 - 2*m2*l1*l2*p[1]*p[2]*cos(θ[1]-θ[2]) ) /
        (2*m2*l1^2*l2^2*(m1+m2*sin(θ[1]-θ[2])^2)) -
        (m1+m2)*g*l1*cos(θ[1]) - m2*g*l2*cos(θ[2])
end

l1 = 1. #length of pendulum1
l2 = 2. #length of pendulum2
m1 = 1. #mass of pendulum1
m2 = 1. #mass of pendulum2
g = 9.81 #gravity

params = [l1, l2, m1, m2, g]
q0 = [1.0,1.0]
p0 = [1.0,1.0]
times = (0.,25.)
prob = HamiltonianProblem(H, q0, p0, times, params)
sol2 = solve(prob, SofSpa10(), dt = .05)

plot(sol2, vars=1, xlim=(0,20), label="Momentum1")
plot!(sol2, vars=2, xlim=(0,20), label="Momentum2")
plot!(sol2, vars=3, xlim=(0,20), label="theta1")
plot!(sol2, vars=4, xlim=(0,20), label="theta2")
#==========================================================#


#Animation
#==========================================================#
function make_pretty_gif(sol)
    timepoints = sol.t

    x1 = l1*sin.(sol[3,:])
    y1 = -l1*cos.(sol[3,:])
    x2 = x1 + l2*sin.(sol[4,:])
    y2 = y1 - l2*cos.(sol[4,:])

    axis_lim = (l1+l2)*1.2

    anim = Animation()
    for i =1:length(timepoints)
        str = string("Time = ", round(timepoints[i],digits=1), " sec")
        plot([0,x1[i]], [0,y1[i]], size=(400,300), xlim=(-axis_lim,axis_lim), ylim=(-axis_lim,1), markersize = 10, markershape = :circle,label ="",axis = [])
        plot!([x1[i],x2[i]], [y1[i],y2[i]], markersize = 10, markershape = :circle,label ="",title = str, title_location = :left)

        if i > 8 #rainbow trail
            plot!([x2[i-2:i]],   [y2[i-2:i]],  alpha = 0.15, linewidth = 2, color = :red, label=nothing)
            plot!([x2[i-3:i-2]], [y2[i-3:i-2]],alpha = 0.15, linewidth = 2, color = :orange, label=nothing)
            plot!([x2[i-4:i-3]], [y2[i-4:i-3]],alpha = 0.15, linewidth = 2, color = :yellow, label=nothing)
            plot!([x2[i-6:i-4]], [y2[i-6:i-4]],alpha = 0.15, linewidth = 2, color = :green, label=nothing)
            plot!([x2[i-7:i-6]], [y2[i-7:i-6]],alpha = 0.15, linewidth = 2, color = :blue, label=nothing)
            plot!([x2[i-8:i-7]], [y2[i-8:i-7]],alpha = 0.15, linewidth = 2, color = :purple, label=nothing)
        end
        frame(anim)
    end
    gif(anim, fps = 30)
end

make_pretty_gif(sol1)
make_pretty_gif(sol2)
