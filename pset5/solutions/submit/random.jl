function randomUnif(seed, n)
	d1 = 2147483647.0
	d2 = 2147483711.0
	dmultix = 16807.0

	r = zeros(n)

	for i in 1:n
		seed = dmultix*seed % d1
		r[i] = seed / d2

	end

	if (n ==1)
		return r[1], seed
	else
		return r, seed
	end
end

