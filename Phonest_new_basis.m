%question same for both banks 4-states
% specifically working with XX
clear all;
tic;

%4-dim orthonormal vector
trials = 2;

alpha_homo_cell = cell(1,trials);
alpha2_homo_cell = cell(1,trials);
homodyne_cell = cell(1,trials);

for t = 0: trials
    for c = 1:5
        m(c) = (0.2*c);
        m2(c) = (0.2*t) + m(c);

        alpha = m(c);
        alpha2 = m2(c);

        split = (alpha + alpha2)/2;
        
        if alpha2 < 1
            n = 4;
        elseif alpha2 < 1.3
            n = 6;
        else 
            n=8;
        end

        %n=4;
        fprintf('n = %d\n', n);
          % Truncation level
        tolerance = 1e-10;     % Tolerance for zero norm


        mpi0 =  zeros(n);
        for i=1:n
            for j=1:n
                mpi0(i,j) = pi0(i-1,j-1,split);
            end
        end
        
        mpi1 =  zeros(n);
        for i=1:n
            for j=1:n
                mpi1(i,j) = pi1(i-1,j-1,split);
            end
        end
        
        mpi2 =  zeros(n);
        for i=1:n
            for j=1:n
                mpi2(i,j) = pi2(i-1,j-1,split);
               % mpi2(i,j) = (-1)^(i+j-2)*pi0(i-1,j-1,split);
            end
        end
        
        
        mpi3 =  zeros(n);
        for i=1:n
            for j=1:n
                mpi3(i,j) = pi3(i-1,j-1,split);
                %mpi3(i,j) = (-1)^(i+j-2)*pi1(i-1,j-1,split);
            end
        end

        %mpi0+mpi2 = identity
        %mpi1+mpi3 = identity

        %mpi00 = kron(mpi0,mpi0);
        %mpi11 = kron(mpi1,mpi1);
        %mpi22 = kron(mpi2,mpi2);
        %mpi33 = kron(mpi3,mpi3);


        coherent_params = [alpha, 1i * alpha, -alpha, -1i * alpha,alpha2,1i * alpha2, -alpha2,-1i * alpha2 ];


        % --- Generate initial state vectors ---
      
        initial_vectors = cell(1, length(coherent_params));
        for k = 1:length(coherent_params)
            initial_vectors{k} = get_fock_vector(coherent_params(k), n, tolerance);
        end 

        %print()


        fprintf('Generated %d initial vectors in Fock basis (dim=%d).\n', length(initial_vectors), n);


        % --- Display individual vectors ---
        for i = 1:8
            fprintf('Initial coherent |alpha_%d> = \n', i-1);
            disp(initial_vectors{i});
        end


        fprintf('\nStarting Numerical Gram-Schmidt process...\n');
        psi_vectors = gram_schmidt_numeric(initial_vectors, tolerance);
        fprintf('Gram-Schmidt process completed.\n');
   
        
        % --- Display individual vectors ---
        for i = 1:8
            fprintf('|psi_%d> = \n', i-1);
            disp(psi_vectors{i});
        end

        num_basis = length(psi_vectors);
        num_coherent = length(initial_vectors);

        Psi = cell2mat(psi_vectors);  % N_max x num_basis
        A = cell2mat(initial_vectors); % N_max x num_coherent

        B = Psi' * A;  % B(k, j) = <psi_k | alpha_j>

        % Optional: verify reconstruction (only if num_basis == num_coherent)
        if num_basis == num_coherent
            A_reconstructed = Psi * B;  % Should be close to A
            error = norm(A - A_reconstructed);  % Should be small.
            disp("Reconstruction error (should be ~0):")
            disp(error)
        else
            disp("Cannot fully reconstruct - basis is incomplete")
        end
        
        alpha_in_psi_basis = cell(1, num_coherent);
        
        for j = 1:num_coherent
            alpha_in_psi_basis{j} = B(:, j);  % 8x1 vector: coefficients of alpha_j in psi_k basis
        end
        

        fprintf('\nOrthonormal basis vectors (|psi_k>) in Fock basis:\n');
        for k = 1:length(alpha_in_psi_basis)
            norm_val = norm(alpha_in_psi_basis{k});
            fprintf('|alpha_%d>: shape=[%d], norm=%.4f\n', k-1, length(alpha_in_psi_basis{k}), norm_val);
            fprintf('%+.4f\n', alpha_in_psi_basis{k});
        end
        
        % --- Verify Orthogonality ---
        fprintf('\nVerifying orthogonality <psi_i | psi_j>:\n');
        orthogonality_check = true;
        for i = 1:length(alpha_in_psi_basis)
            for j = 1:length(alpha_in_psi_basis)
                inner_product = alpha_in_psi_basis{i}' * alpha_in_psi_basis{j};
                expected = (i == j);
                is_close = abs(inner_product - expected) < tolerance;
                fprintf('<alpha_%d | alpha_%d> = %.3f (Expected: %d) -> Correct: %d\n', ...
                    i-1, j-1, real(inner_product), expected, is_close);
                if ~is_close
                    orthogonality_check = false;
                end
            end 
        end
        
        if orthogonality_check
            fprintf('\nOrthogonality verification successful (within tolerance).\n');
        else
            fprintf('\nWARNING: Orthogonality verification failed for some pairs.\n');
        end
        
        % --- Display individual vectors ---
        for i = 1:8
            fprintf('|alpha_%d> = \n', i-1);
            disp(alpha_in_psi_basis{i});
        end 
        
        %{
        %Density Operator
        dc0 = B(:,1) * ctranspose(B(:,1));
        dc1 = B(:,2) * ctranspose(B(:,2));
        dc2 = B(:,3) * ctranspose(B(:,3));
        dc3 = B(:,4) * ctranspose(B(:,4));
        dc4 = B(:,5) * ctranspose(B(:,5));
        dc5 = B(:,6) * ctranspose(B(:,6));
        dc6 = B(:,7) * ctranspose(B(:,7));
        dc7 = B(:,8) * ctranspose(B(:,8));
        %}
        full_state_0 = Psi * B(:,1);  % This is the |alpha_0> state in Fock basis

        dc0 = full_state_0 * full_state_0';
        dc1 = (Psi * B(:,2)) * (Psi * B(:,2))';
        dc2 = (Psi * B(:,3)) * (Psi * B(:,3))';
        dc3 = (Psi * B(:,4)) * (Psi * B(:,4))';
        dc4 = (Psi * B(:,5)) * (Psi * B(:,5))';
        dc5 = (Psi * B(:,6)) * (Psi * B(:,6))';
        dc6 = (Psi * B(:,7)) * (Psi * B(:,7))';
        dc7 = (Psi * B(:,8)) * (Psi * B(:,8))';

        eta5(c) =1/4*(trace(dc0*mpi0) + trace(dc4*mpi1) + trace(dc2*mpi2) + trace(dc6*mpi3)); %(so system: AB)


    end

    alpha_homo_cell{t+1} = m;
    alpha2_homo_cell{t+1} = m2;
    homodyne_cell{t+1} = eta5;
    save('alpha_homo_cell.mat', 'alpha_homo_cell');
    save('alpha2_homo_cell.mat', 'alpha2_homo_cell');
    save('homodyne_cell.mat', 'homodyne_cell');
 % alpha2_all{t+1} = m2;
  %eta_all{t+1} = eta1;  % eta_all{run} stores eta1 for this run
end  

    figure;
    hold on;
    
    for t = 1:length(alpha_homo_cell)
        alpha1 = alpha_homo_cell{t};
        alpha2 = alpha2_homo_cell{t};
        eta = homodyne_cell{t};  % eta now varies with alpha2
    
        if length(alpha1) == length(alpha2) && length(alpha2) == length(eta)
            plot3(alpha2, alpha1, eta, 'LineWidth', 1.5);  % z varies with alpha2
        else
            warning('Length mismatch at trial %d: skipping plot.', t);
        end
    end
    
    hold off;
    grid on;
    xlabel('Alpha 2');
    ylabel('Alpha 1');
    zlabel('Probability (eta1)');
    title('3D Line Plot of eta1 varying with alpha2');
    view(3);

        
function f= pi0(n,m,split)
    % n and m \in positive integers and zero
    syms x
    g = hermiteH(n,x)*hermiteH(m,x)*exp(-x^2);
    f1 = int(g,x,0,split);
    f2 = 1.0/sqrt(2^(n+m)*factorial(n)*factorial(m)*pi);
    f = double(f1*f2);
end

function f= pi1(n,m,split)
    % n and m \in positive integers and zero
    syms x
    g = hermiteH(n,x)*hermiteH(m,x)*exp(-x^2);
    f1 = int(g,x,split,inf);
    f2 = 1.0/sqrt(2^(n+m)*factorial(n)*factorial(m)*pi);
    f = double(f1*f2);
end

function f= pi2(n,m,split)
    % n and m \in positive integers and zero
    syms x
    g = hermiteH(n,x)*hermiteH(m,x)*exp(-x^2);
    f1 = int(g,x,-split,0);
    f2 = 1.0/sqrt(2^(n+m)*factorial(n)*factorial(m)*pi);
    f = double(f1*f2);
end

function f= pi3(n,m,split)
    % n and m \in positive integers and zero
    syms x
    g = hermiteH(n,x)*hermiteH(m,x)*exp(-x^2);
    f1 = int(g,x,-inf,-split);
    f2 = 1.0/sqrt(2^(n+m)*factorial(n)*factorial(m)*pi);
    f = double(f1*f2);
end

function H = hermiteH(n, x)
% Physicist's Hermite polynomial H_n(x)
    if n == 0
        H = ones(size(x));
    elseif n == 1
        H = 2 * x;
    else
        Hn_2 = ones(size(x));
        Hn_1 = 2 * x;
        for k = 2:n
            H = 2 * x .* Hn_1 - 2 * (k - 1) * Hn_2;
            Hn_2 = Hn_1;
            Hn_1 = H;
        end
    end
end

function vec = get_fock_vector(alpha, N_max, tolerance)
    abs_alpha_sq = abs(alpha)^2;
    prefactor = exp(-abs_alpha_sq / 2.0);

    vec = zeros(N_max, 1);
    n = 0:(N_max - 1);
    sqrt_n_fact = sqrt(factorial(n));

    if alpha == 0
        vec(1) = 1.0;
        return;
    end

    alpha_pow_n = alpha.^n;
    vec = prefactor * (alpha_pow_n ./ sqrt_n_fact).';
    
    vec(abs(vec) < tolerance) = 0.0;
end


% --- Gram-Schmidt Orthonormalization ---
function basis = gram_schmidt_numeric(vectors, tol)
    num_vectors = length(vectors);
    basis = cell(1, num_vectors);

    for i = 1:num_vectors
        fprintf('Processing vector %d...\n', i);
        vi = vectors{i};

        for j = 1:i-1
            proj_coeff = basis{j}' * vi;
            vi = vi - proj_coeff * basis{j};
        end

        norm_vi = norm(vi);
        if norm_vi < tol
            fprintf('Warning: Vector %d became (close to) zero. Norm = %e.\n', i, norm_vi);
            basis{i} = zeros(size(vi));
        else
            basis{i} = vi / norm_vi;
        end
        fprintf('Finished vector %d. Current basis size: %d\n', i, i);
    end
end