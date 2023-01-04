
# Standardized version of the location/time variables.
Yz = copy(Y)
for j in 1:3
    Yz[:, j] .-= mean(Yz[:, j])
    Yz[:, j] ./= std(Yz[:, j])
end


function slicemeans(X, nslice)
    n, p = size(X)
    mn = zeros(p, nslice)
    ss = Int(ceil(n / nslice)) # Slice size
    wt = zeros(nslice)
    jj = 0
    for i in 1:nslice
        j2 = min(jj+ss, n)
        mn[:, i] = mean(X[jj+1:j2, :], dims=1)[:]
        wt[i] = j2 - jj
        jj += ss
    end
    wt ./= sum(wt)

    c = zeros(p, p)
    for i in eachindex(wt)
        c .+= wt[i] * mn[:, i] * mn[:,i]'
    end

    return c
end

function mfsir(Y, X; k=5, nslice=20, ndir=3)

    # Center the data
    n, p = size(X)
    xmn = mean(X, dims=1)[:]
    for j in 1:p
        X[:, j] .-= xmn[j]
    end

    # Get the marginal covariance
    mcov = Symmetric(cov(X))

    # Zero out the low-order eigenvectors
    # of the marginal covariance
    eg = eigen(mcov)
    Pf = eg.vectors[:, end-k+1:end]
    mcovp = Pf' * mcov * Pf
    mcovp = Pf * mcovp * Pf'
    mcovp = Symmetric(mcovp)

    # Prepare to take the inverse square root
    # of mcovp.
    eg = eigen(mcovp)
    ii = findall(eg.values .> 1e-5)
    w = sqrt.(eg.values[ii])
    wi = 1 ./ w
    v = eg.vectors[:, ii]

    y = Y * randn(size(Y, 2))
    ii = sortperm(y)
    y = y[ii]
    X = X[ii, :]

    scov = slicemeans(X, nslice)

    # Transform the slice covariance to
    # decorrelated coordinates.
    scovs = v * Diagonal(wi) * v' * scov * v * Diagonal(wi) * v'

    # Get the basis vectors in the decorrelated coordinates,
    # then transform back to the original coordinates.
    eg = eigen(Symmetric(scovs))
    ii = sortperm(eg.values, rev=true)
    w = eg.values[ii]
    u = eg.vectors[:, ii]
    b = v * Diagonal(wi) * v' * u[:, 1:ndir]

    return b
end





b = mfsir(Y, temp)

S = tempc * b
for j in 1:3
    for k in 1:3
        m = loess(S[:, j], Y[:, k])
        us = range(extrema(S[:, j])..., 100)
        ss = predict(m, us)
        plt = lineplot(us, ss)
        println(plt)
    end
end
