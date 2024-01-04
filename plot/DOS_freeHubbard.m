% numeric estimate the DOS for free electron gas with only NN hopping

L = 10000; % lattice size
N = L * L;
k_set = -pi:2*pi/L:pi-2*pi/L;
doping = 1/8;
dE = 0.0001;

% calculate chemical potential
energy_map = zeros(L);
assert(numel(k_set) == L)
for i = 1:L
    kx = k_set(i);
    energy_map(i, :) = get_energy(kx,k_set);
end
energy_list = reshape(energy_map,[],1);
energy_list = sort(energy_list);
mu = energy_list((1-doping) * 1/2 * N);

fprintf("mu = %.6f\n", mu);

energy_select = energy_list(energy_list > mu - dE & energy_list < mu+dE);

num_of_state = numel(energy_select);
dos = num_of_state/N /(2*dE) * 2; % *2 for spin dof
fprintf("DOS = %.6f\n", dos);


function e = get_energy(kx,ky)
e = -2 * (cos(kx) + cos(ky));
end