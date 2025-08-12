%question same for both banks 4-states
% specifically working with XX
clear all;
tic;

%4-dim orthonormal vector
trials = 2;

alphaxx_cell = cell(1,trials);
alpha2xx_cell = cell(1,trials);
eta1xx_cell = cell(1,trials);

for t = 0: trials
    for c = 1:5
        m(c) = (0.2*c);
        m2(c) = (0.2*t) + m(c);

        alpha = m(c);
        alpha2 = m2(c);

        % Truncation level
        %N_max= 8; In case you want to set n with a static value 
        alpha_abs = abs(alpha2);
        N_max = ceil(alpha_abs^2 + 6 * sqrt(alpha_abs^2));

        % Tolerance for zero norm
        tolerance = 1e-10;     


        coherent_params = [alpha, 1i * alpha, -alpha, -1i * alpha,alpha2,1i * alpha2, -alpha2,-1i * alpha2 ];


        % --- Generate initial state vectors ---
      
        initial_vectors = cell(1, length(coherent_params));
        for k = 1:length(coherent_params)
            initial_vectors{k} = get_fock_vector(coherent_params(k), N_max, tolerance);
        end 


        fprintf('Generated %d initial vectors in Fock basis (dim=%d).\n', length(initial_vectors), N_max);


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
        

    
        %Density Operator
        dc0 = B(:,1) * ctranspose(B(:,1));
        dc1 = B(:,2) * ctranspose(B(:,2));
        dc2 = B(:,3) * ctranspose(B(:,3));
        dc3 = B(:,4) * ctranspose(B(:,4));
        dc4 = B(:,5) * ctranspose(B(:,5));
        dc5 = B(:,6) * ctranspose(B(:,6));
        dc6 = B(:,7) * ctranspose(B(:,7));
        dc7 = B(:,8) * ctranspose(B(:,8));
     
        %SDP    
        cvx_begin sdp quiet

            variable X01(8,8) complex semidefinite;
            variable X12(8,8) complex semidefinite;
            variable X23(8,8) complex semidefinite;
            variable X30(8,8) complex semidefinite;

            maximize real(0.25*trace(dc0*X01+dc4*X12 +dc2*X23 +dc6*X30)) %objective function
                X01+X12+X23+X30 == eye(8);
        cvx_end
        cvx_status
        eta1(c) = real(0.25*trace(dc0*X01+dc4*X12 +dc2*X23 +dc6*X30));

        %eta1(c) = 0.5*(eta1(c)+1); %Then normalization that's in the equation

    end
    alphaxx_cell{t+1} = m;
    alpha2xx_cell{t+1} = m2;
    eta1xx_cell{t+1} = eta1;
    save('alphaxx_cell.mat', 'alphaxx_cell');
    save('alpha2xx_cell.mat', 'alpha2xx_cell');
    save('eta1xx_cell.mat', 'eta1xx_cell');
end  


    figure;
    hold on;
    
    for t = 1:length(alphaxx_cell)
        alpha1 = alphaxx_cell{t};
        alpha2 = alpha2xx_cell{t};
        eta = eta1xx_cell{t};  % eta now varies with alpha2
    
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