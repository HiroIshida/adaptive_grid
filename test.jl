using Interpolations
include("octree.jl")

function test_2dim(tol)
    sigma = 7
    R(a) = [cos(a) -sin(a);
            sin(a) cos(a)]
    function sdf(x, a, b)
        x = R(-a)*x
        d = abs.(x) - b
        tmp = [max(d[1], 0.0), max(d[2], 0.0)]
        return norm(tmp) + min(max(d[1], d[2]), 0.0)
    end

    a = pi/10
    b = [60, 30]
    f(x) = 0.5*(1 + erf(sdf(x, a, b)/sqrt(2*sigma^2)))
    n_grid = 20
    predicate(node::Node) = pred_standard(node, f, tol, n_grid, Linear())

    tree = Tree([-100, -100], [100, 100], f)
    auto_split!(tree, predicate)
    vertex_reduction!(tree)
    #show(tree)
    
    isValid = true
    for i in 1:1000
        myrn() = rand()*200 - 100
        q = [myrn(), myrn()]
        error = abs(evaluate(tree, q) - f(q))
        if error > tol*100
            isValid *= false
        end
    end
    return isValid
end

function test_3dim(tol)
    sigma = 7
    function sdf(x, b)
        d = abs.(x) - b
        tmp = [max(d[1], 0.0), max(d[2], 0.0), max(d[3], 0.0)]
        return norm(tmp) + min(max(d[1], d[2], d[3]), 0.0)
    end

    b = [60, 30, 30]
    f(x) = 0.5*(1 + erf(sdf(x, b)/sqrt(2*sigma^2)))
    n_grid = 1
    predicate(node::Node) = pred_standard(node, f, tol, n_grid, Linear())


    tree = Tree([-100, -100, -100], [100, 100, 100], f)
    auto_split!(tree, predicate)
    vertex_reduction!(tree)
    
    isValid = true
    for i in 1:1000
        myrn() = rand()*200 - 100
        q = [myrn(), myrn(), myrn()]
        error = abs(evaluate(tree, q) - f(q))
        if error > tol*10
            isValid *= false
        end
    end
    return isValid
end

tol = 0.1
test_2dim(tol)
test_3dim(tol)
