%% LaTeX Expression Generator for Versal Deformation
% This script generates LaTeX expressions for all matrices and formulas
% in the versal deformation construction algorithm
% Based on ALP 2026 and SIAM 2003 (Garcia-Planas & Mailybaev)
%
% UPDATED: Numerical stabilization via SVD-based orthonormalization
%          (preserves T*Vs = St and makes K = I)

%% Initial setup
format compact
digits(3); % 3-digit precision
syms a0 a1 a2 real; % Symbolic coefficients

% Numerical values for calculation (quadratic polynomial example)
a0v = 1; a1v = 2; a2v = 3;

%% Step 0: Build companion matrices
[Ao, Bo] = step0_build_Ao_Bo(a0v, a1v, a2v);
latex_Ao = latex(sym(Ao));
latex_Bo = latex(sym(Bo));

%% Step 1: Build matrix T (tangent space representation)
T = step1_build_T(a0v, a1v, a2v);
latex_T = latex(vpa(sym(T)));

%% Step 2: Build constraint matrix T1 (optimal subspace conditions)
T1 = step2_build_T1(a0v, a1v, a2v);
latex_T1 = latex(vpal(sym(T1),3));

%% Step 3: Kernel of T1
kerT1 = step3_kernel_T1(T1);
latex_kerT1 = latex(vpal(sym(kerT1),3));

%% Step 4: Build orthonormal basis of ker(T1)
Z = orth(kerT1);               % orthonormal basis of ker(T1)
A = T * Z;                     % apply T

%% Step 5: SVD of A to construct St and Vs (ensures K = I)
[U, Sigma, V] = svd(A, 'econ');
tol = max(size(A)) * eps(max(diag(Sigma)));
keep = diag(Sigma) > tol;
St = U;                         % orthonormal basis of S_t
Sigma_inv = diag(1 ./ diag(Sigma(keep,keep)));
Vs = Z * V(:,keep) * Sigma_inv; % such that T * Vs = St

% Now St' * St = I, and T * Vs = St, so K = St' * (T * Vs) = I
latex_St = latex(vpal(sym(St),3));
latex_Vs = latex(vpal(sym(Vs),3));

%% Step 6: Orthogonal complement of St (within the full pencil space)
Ncomp = null(St');   % orthonormal basis of St^⊥
latex_Ncomp = latex(vpal(sym(Ncomp),3));

%% Step 7: Orthogonal complement of Vs (optional)
Vscomp = null(Vs');  % orthonormal basis of Vs^⊥
latex_Vscomp = latex(vpal(sym(Vscomp),3));

%% Step 8: Basis of generalized Sylvester space GS
GS = step7_generalized_Sylvester();
latex_GS = latex(sym(GS));

%% Step 9: Matrices L and K for recurrence relations
L = GS' * Ncomp;               % 3x3 matrix
K = St' * (T * Vs);            % Should be exactly identity (up to round-off)
latex_L = latex(vpal(sym(L),3));
latex_K = latex(vpal(sym(K),3));

%% Step 10: First and second-order derivatives
[phi1, mu1, Bphi, Bmu] = step9_derivatives(L, K, St, Ncomp, GS, Vs, Ao, Bo);

% Matrix A for first-order derivatives of phi
latex_A_phi = latex(vpal(sym(phi1),3));

% Matrix A^mu for first-order derivatives of mu
latex_A_mu = latex(vpal(sym(mu1),3));

% Matrices B for second-order derivatives of phi
B0 = squeeze(Bphi(1,:,:));
B1 = squeeze(Bphi(2,:,:));
B2 = squeeze(Bphi(3,:,:));
latex_B0 = latex(vpal(sym(B0),3));
latex_B1 = latex(vpal(sym(B1),3));
latex_B2 = latex(vpal(sym(B2),3));

% Matrices B^mu for second-order derivatives of mu
Bmu0 = squeeze(Bmu(1,:,:));
Bmu1 = squeeze(Bmu(2,:,:));
Bmu2 = squeeze(Bmu(3,:,:));
Bmu3 = squeeze(Bmu(4,:,:));
Bmu4 = squeeze(Bmu(5,:,:));
latex_Bmu0 = latex(vpal(sym(Bmu0),3));
latex_Bmu1 = latex(vpal(sym(Bmu1),3));
latex_Bmu2 = latex(vpal(sym(Bmu2),3));
latex_Bmu3 = latex(vpal(sym(Bmu3),3));
latex_Bmu4 = latex(vpal(sym(Bmu4),3));

% Generate LaTeX formulas
latex_formulas = generate_versal_formulas(phi1, Bphi);
latex_mu_formulas = generate_mu_formulas(mu1, Bmu);
latex_pencil_formulas = generate_pencil_formulas(Vs);

%% Generate LaTeX output file (normalized bases)
fid = fopen('versal_deformation_output.tex', 'w');

fprintf(fid, '\\documentclass{article}\n');
fprintf(fid, '\\usepackage{amsmath, amssymb}\n');
fprintf(fid, '\\usepackage{geometry}\n');
fprintf(fid, '\\geometry{a4paper, margin=1in}\n');
fprintf(fid, '\\title{Versal Deformation Results (Normalized Bases)}\n');
fprintf(fid, '\\begin{document}\n');
fprintf(fid, '\\maketitle\n');

% Step 0: Companion matrices
fprintf(fid, '\\section*{Step 0: Companion Matrices}\n');
fprintf(fid, 'Matrix $A_o$: \n\\[\n%s\n\\]\n', latex_Ao);
fprintf(fid, 'Matrix $B_o$: \n\\[\n%s\n\\]\n', latex_Bo);

% Step 1: Matrix T
fprintf(fid, '\\section*{Step 1: Matrix T}\n');
fprintf(fid, '\\[\nT = %s\n\\]\n', latex_T);

% Step 2: Matrix T1
fprintf(fid, '\\section*{Step 2: Matrix T1}\n');
fprintf(fid, '\\[\nT_1 = %s\n\\]\n', latex_T1);

% Step 3: Kernel of T1
fprintf(fid, '\\section*{Step 3: Basis of $\\mathrm{Ker} T_1$}\n');
fprintf(fid, '\\[\n\\mathrm{Ker} T_1 = %s\n\\]\n', latex_kerT1);

% Step 4: Basis of St (normalized)
fprintf(fid, '\\section*{Step 4: Basis of $S_t$ (normalized)}\n');
fprintf(fid, '\\[\nS_t = %s\n\\]\n', latex_St);

% Step 5: Orthogonal complement of St
fprintf(fid, '\\section*{Step 5: Basis of $S_t^\\perp$}\n');
fprintf(fid, '\\[\nS_t^\\perp = %s\n\\]\n', latex_Ncomp);

% Step 6: Basis of Vs (normalized)
fprintf(fid, '\\section*{Step 6: Basis of $\\mathcal{V}_s$ (normalized)}\n');
fprintf(fid, '\\[\n\\mathcal{V}_s = %s\n\\]\n', latex_Vs);

% Step 6a: Orthogonal complement of Vs
fprintf(fid, '\\section*{Step 6a: Orthogonal Complement of $\\mathcal{V}_s$}\n');
fprintf(fid, '\\[\n\\mathcal{V}_s^\\perp = %s\n\\]\n', latex_Vscomp);

% Step 7: Basis of GS
fprintf(fid, '\\section*{Step 7: Basis of $\\mathcal{GS}$}\n');
fprintf(fid, '\\[\n\\mathcal{GS} = %s\n\\]\n', latex_GS);

% Step 8: Matrices L and K
fprintf(fid, '\\section*{Step 8: Matrices L and K}\n');
fprintf(fid, 'Matrix $L$: \n\\[\nL = %s\n\\]\n', latex_L);
fprintf(fid, 'Matrix $K$: \n\\[\nK = %s\n\\]\n', latex_K);

% Step 9: First-order derivatives
fprintf(fid, '\\section*{Step 9: First-Order Derivatives}\n');
fprintf(fid, 'Matrix $A$ ($\\varphi$ derivatives):\n\\[\nA = %s\n\\]\n', latex_A_phi);
fprintf(fid, 'Matrix $A^\\mu$ ($\\mu$ derivatives):\n\\[\nA^\\mu = %s\n\\]\n', latex_A_mu);

% Step 9: Second-order derivatives of phi
fprintf(fid, '\\section*{Step 9: Second-Order Derivatives of $\\varphi$}\n');
fprintf(fid, 'Matrix $B^0$ (for $\\varphi_0$):\n\\[\nB^0 = %s\n\\]\n', latex_B0);
fprintf(fid, 'Matrix $B^1$ (for $\\varphi_1$):\n\\[\nB^1 = %s\n\\]\n', latex_B1);
fprintf(fid, 'Matrix $B^2$ (for $\\varphi_2$):\n\\[\nB^2 = %s\n\\]\n', latex_B2);

% Step 9: Second-order derivatives of mu
fprintf(fid, '\\section*{Step 9: Second-Order Derivatives of $\\mu$}\n');
fprintf(fid, 'Matrix $B^{\\mu 0}$ (for $\\mu_1$):\n\\[\nB^{\\mu 0} = %s\n\\]\n', latex_Bmu0);
fprintf(fid, 'Matrix $B^{\\mu 1}$ (for $\\mu_2$):\n\\[\nB^{\\mu 1} = %s\n\\]\n', latex_Bmu1);
fprintf(fid, 'Matrix $B^{\\mu 2}$ (for $\\mu_3$):\n\\[\nB^{\\mu 2} = %s\n\\]\n', latex_Bmu2);
fprintf(fid, 'Matrix $B^{\\mu 3}$ (for $\\mu_4$):\n\\[\nB^{\\mu 3} = %s\n\\]\n', latex_Bmu3);
fprintf(fid, 'Matrix $B^{\\mu 4}$ (for $\\mu_5$):\n\\[\nB^{\\mu 4} = %s\n\\]\n', latex_Bmu4);

% Versal deformation formulas
fprintf(fid, '\\section*{Versal Deformation Formulas}\n');
fprintf(fid, '%s\n', latex_formulas);
fprintf(fid, '\\section*{Versal Deformation Formulas for mu}\n');
fprintf(fid, '%s\n', latex_mu_formulas);
fprintf(fid, '\\section*{Pencil Representation of Basis Vectors}\n');
fprintf(fid, '%s\n', latex_pencil_formulas);

fprintf(fid, '\\end{document}');
fclose(fid);

disp('LaTeX file generated: versal_deformation_output.tex');

%% Generate comparison LaTeX (raw vs normalized) for supplementary material
% For completeness, compute raw bases (non-orthonormal) for comparison
St_raw = T * kerT1;
Vs_raw = Vs;   % but Vs is already from SVD; we can compute a raw version:
% Actually Vs_raw from projection method (if needed) - but not essential.
% We'll just show the normalized ones as main result.

fid2 = fopen('bases_comparison_supplement.tex', 'w');
fprintf(fid2, '\\documentclass{article}\n');
fprintf(fid2, '\\usepackage{amsmath, amssymb}\n');
fprintf(fid2, '\\usepackage{geometry}\n');
fprintf(fid2, '\\geometry{a4paper, margin=1in}\n');
fprintf(fid2, '\\title{Comparison of Non-Normalized and Normalized Bases}\n');
fprintf(fid2, '\\begin{document}\n');
fprintf(fid2, '\\maketitle\n');
fprintf(fid2, '\\section*{Non-Normalized Bases (raw)}\n');
fprintf(fid2, '$S_t^{\\text{raw}} = %s$\n\n', latex(vpal(sym(St_raw),3)));
% For Vs_raw we can compute a simple non-orthogonal preimage: Vs_raw_simple = kerT1;
Vs_raw_simple = kerT1;
fprintf(fid2, '$\\mathcal{V}_s^{\\text{raw}} = %s$\n\n', latex(vpal(sym(Vs_raw_simple),3)));
fprintf(fid2, '\\section*{Normalized Bases (orthonormal, $K=I$)}\n');
fprintf(fid2, '$S_t^{\\text{norm}} = %s$\n\n', latex_St);
fprintf(fid2, '$\\mathcal{V}_s^{\\text{norm}} = %s$\n\n', latex_Vs);
fprintf(fid2, '\\end{document}\n');
fclose(fid2);
disp('Comparison LaTeX file (raw vs normalized) generated: bases_comparison_supplement.tex');

%% ========================================================================
%% AUXILIARY FUNCTIONS (same as yours)
%% ========================================================================

function [Ao, Bo] = step0_build_Ao_Bo(a0, a1, a2)
    Ao = [a1, a0; -1, 0];
    Bo = [-a2, 0; 0, -1];
end

function T = step1_build_T(a0, a1, a2)
    Ao = [a1, a0; -1, 0];
    Bo = [-a2, 0; 0, -1];
    p = size(Ao, 1);
    q = size(Ao, 2);
    T = [kron(Ao', eye(p)), -kron(eye(q), Ao);
         kron(Bo', eye(p)), -kron(eye(q), Bo)];
end

function T1 = step2_build_T1(~, ~, ~)
    T1 = [0,0,0,0, 0,1,0,0;    % v21 = 0
          0,0,0,0, 0,0,0,1;    % v22 = 0
          1,0,0,1, -1,0,0,-1]; % trace condition
end

function ker = step3_kernel_T1(T1)
    ker = null(T1);
end

function GS = step7_generalized_Sylvester()
    GS = [1, 0, 0, 0, 0, 0, 0, 0;
          0, 0, 1, 0, 0, 0, 0, 0;
          0, 0, 0, 0, -1, 0, 0, 0]';
end

function [phi1, mu1, Bphi, Bmu] = step9_derivatives(L, K, St, Ncomp, GS, Vs, Ao, Bo)
    len = 8;
    phi1 = zeros(size(GS,2), len);
    mu1  = zeros(size(Vs,2), len);
    
    % First-order derivatives
    for k = 1:len
        ek = zeros(len,1); ek(k)=1;
        rhs = Ncomp' * ek;
        phi1(:,k) = L \ rhs;
        temp = ek - GS * phi1(:,k);
        rhs2 = St' * temp;
        mu1(:,k) = K \ rhs2;
    end
    
    % Second-order derivatives (using pinv for stability, but K should be I)
    c = mat2cell(GS, size(GS,1), ones(1, size(GS,2)));
    t = mat2cell(St, size(St,1), ones(1, size(St,2)));
    r = mat2cell(Vs, size(Vs,1), ones(1, size(Vs,2)));
    n = mat2cell(Ncomp, size(Ncomp,1), ones(1, size(Ncomp,2)));
    
    Z = L;
    A = pinv(Z, 1e-10) * Ncomp';
    mu = cell(1,8);
    for i = 1:8
        fb = zeros(8,1); fb(i)=1;
        expr = fb - A(1,i)*c{1} - A(2,i)*c{2} - A(3,i)*c{3};
        b_vals = zeros(1,5);
        for kk=1:5
            b_vals(kk) = eprodu1(expr, t{kk});
        end
        mu{i} = pinv(K,1e-10) * b_vals';
    end
    gmu = [mu{1}, mu{2}, mu{3}, mu{4}, mu{5}, mu{6}, mu{7}, mu{8}];
    
    alfaij = cell(8,8);
    for i = 1:8
        for j = 1:8
            X_i = A(1,i)*c{1} + A(2,i)*c{2} + A(3,i)*c{3};
            Y_j = gmu(1,j)*r{1} + gmu(2,j)*r{2} + gmu(3,j)*r{3} + gmu(4,j)*r{4} + gmu(5,j)*r{5};
            Z_i = zeros(8,1); Z_i(i)=1;
            alfaij{i,j} = -2 * galfa_three_correcto(X_i, Y_j, Z_i, Ao, Bo);
        end
    end
    
    S2 = zeros(8,1);
    for i=1:8, for j=1:8, S2 = S2 + alfaij{i,j}; end; end
    
    MM = zeros(3,64);
    for k=1:3
        for i=1:8
            for j=1:8
                idx = 8*(i-1)+j;
                MM(k,idx) = eprodu1(S2, n{k});
            end
        end
    end
    F = pinv(Z,1e-10);
    Dphi = F * MM;
    Bphi = zeros(3,8,8);
    for k=1:3, Bphi(k,:,:) = reshape(Dphi(k,:),8,8); end
    
    d2mu = cell(1,8);
    for i=1:8
        fb = zeros(8,1); fb(i)=1;
        expr = fb - squeeze(Bphi(1,i,i))*c{1} - squeeze(Bphi(2,i,i))*c{2} - squeeze(Bphi(3,i,i))*c{3};
        b_vals = zeros(1,5);
        for kk=1:5, b_vals(kk) = eprodu1(expr, t{kk}); end
        d2mu{i} = pinv(K,1e-10) * b_vals';
    end
    gmu2 = [d2mu{1}, d2mu{2}, d2mu{3}, d2mu{4}, d2mu{5}, d2mu{6}, d2mu{7}, d2mu{8}];
    
    MMM = zeros(5,64);
    for k=1:5
        for i=1:8
            for j=1:8
                idx = 8*(i-1)+j;
                expr = S2 - squeeze(Bphi(1,i,j))*c{1} - squeeze(Bphi(2,i,j))*c{2} - squeeze(Bphi(3,i,j))*c{3};
                MMM(k,idx) = eprodu1(expr, t{k});
            end
        end
    end
    Fmu = pinv(K,1e-10);
    Dphi2 = Fmu * MMM;
    Bmu = zeros(5,8,8);
    for k=1:5, Bmu(k,:,:) = reshape(Dphi2(k,:),8,8); end
end

function val = galfa_three_correcto(X, Y, Z, Ao, Bo)
    m=2; n=2;
    A_X = reshape(X(1:m*n), m, n);
    B_X = reshape(X(m*n+1:end), m, n);
    U = reshape(Y(1:m*m), m, m);
    V = reshape(Y(m*m+1:end), n, n);
    X_Z = reshape(Z(1:m*n), m, n);
    Y_Z = reshape(Z(m*n+1:end), m, n);
    DeltaA = A_X * V - U * X_Z;
    DeltaB = B_X * V - U * Y_Z;
    val = [DeltaA(:); DeltaB(:)];
end

function val = eprodu1(a, b)
    val = a(:)' * b(:);
end

function B = vpal(A, ndigits, tolerance)
    if nargin < 2, ndigits = 3; end
    if nargin < 3, tolerance = 1e-10; end
    A = vpa(A);
    A(abs(A) < tolerance) = 0;
    B = round(A, ndigits);
end

function latex_output = generate_versal_formulas(phi1, Bphi)
    phi1 = round(phi1, 4);
    Bphi = round(Bphi, 4);
    quad_coeffs = zeros(1,3);
    for k=1:3
        Bk = squeeze(Bphi(k,:,:));
        quad_coeffs(k) = sum(Bk(:)) / 128;
    end
    quad_coeffs = round(quad_coeffs, 3);
    phi_eqs = cell(1,3);
    for k=0:2
        lin_coeffs = phi1(k+1,:);
        lin_terms = {};
        for i=1:8
            coeff = lin_coeffs(i);
            if abs(coeff) < 1e-5, continue; end
            if isempty(lin_terms)
                term = sprintf('%.4f\\epsilon_%d', coeff, i);
            else
                if coeff >= 0, term = sprintf('+ %.4f\\epsilon_%d', coeff, i);
                else term = sprintf('- %.4f\\epsilon_%d', abs(coeff), i); end
            end
            lin_terms{end+1} = term;
        end
        if numel(lin_terms) > 4
            line1 = strjoin(lin_terms(1:min(4,end)), ' ');
            line2 = strjoin(lin_terms(5:end), ' ');
            if ~isempty(line2) && line2(1) ~= '-' && line2(1) ~= '+', line2 = ['+ ' line2]; end
        else
            line1 = strjoin(lin_terms, ' ');
            line2 = '';
        end
        q_coeff = quad_coeffs(k+1);
        if abs(q_coeff) > 0.0005
            quad_sign = '+'; if q_coeff < 0, quad_sign = '-'; end
            quad_term = sprintf('%s %.3f (\\Sigma\\epsilon)^2', quad_sign, abs(q_coeff));
        else
            quad_term = '';
        end
        if isempty(line2)
            eq_str = sprintf('\\varphi_{%d}(\\vec{\\epsilon}) &= %s %s + o(||\\epsilon||^2)', k, line1, quad_term);
        else
            eq_str = sprintf('\\varphi_{%d}(\\vec{\\epsilon}) &= %s \\\\ & \\quad %s %s + o(||\\epsilon||^2)', k, line1, line2, quad_term);
        end
        phi_eqs{k+1} = eq_str;
    end
    latex_output = sprintf(['\\begin{align*}\n%s \\\\\n%s \\\\\n%s\n\\end{align*}\n' ...
        'where $\\Sigma\\epsilon = \\epsilon_{1} + \\epsilon_{2} + \\epsilon_{3} + ' ...
        '\\epsilon_{4} + \\epsilon_{5} + \\epsilon_{6} + \\epsilon_{7} + \\epsilon_{8}$.\n'], ...
        phi_eqs{1}, phi_eqs{2}, phi_eqs{3});
end

function latex_output = generate_mu_formulas(mu1, Bmu)
    mu1 = round(mu1, 3);
    Bmu = round(Bmu, 3);
    quad_coeffs = zeros(1,5);
    for k=1:5
        Bk = squeeze(Bmu(k,:,:));
        quad_coeffs(k) = sum(Bk(:)) / 128;
    end
    quad_coeffs = round(quad_coeffs, 3);
    mu_eqs = cell(1,5);
    for k=1:5
        lin_coeffs = mu1(k,:);
        lin_terms = {};
        for i=1:8
            coeff = lin_coeffs(i);
            if abs(coeff) < 0.0005, continue; end
            if isempty(lin_terms)
                term = sprintf('%.3f\\epsilon_%d', coeff, i);
            else
                if coeff >= 0, term = sprintf('+ %.3f\\epsilon_%d', coeff, i);
                else term = sprintf('- %.3f\\epsilon_%d', abs(coeff), i); end
            end
            lin_terms{end+1} = term;
        end
        linear_str = strjoin(lin_terms, ' ');
        q_coeff = quad_coeffs(k);
        if abs(q_coeff) > 0.0005
            quad_sign = '+'; if q_coeff < 0, quad_sign = '-'; end
            quad_term = sprintf('%s %.3f (\\Sigma\\epsilon)^2', quad_sign, abs(q_coeff));
        else
            quad_term = '';
        end
        if isempty(quad_term)
            eq_str = sprintf('\\mu_{%d}(\\vec{\\epsilon}) = %s + o(||\\epsilon||^2)', k, linear_str);
        else
            eq_str = sprintf('\\mu_{%d}(\\vec{\\epsilon}) = %s \\\\ & %s + o(||\\epsilon||^2)', k, linear_str, quad_term);
        end
        if k < 5, mu_eqs{k} = [eq_str ',']; else mu_eqs{k} = [eq_str '.']; end
    end
    latex_output = sprintf(['\\begin{align*}\n%s\n\\\\\n%s\n\\\\\n%s\n\\\\\n%s\n\\\\\n%s\n\\end{align*}\n'], ...
        mu_eqs{1}, mu_eqs{2}, mu_eqs{3}, mu_eqs{4}, mu_eqs{5});
end

function latex_output = generate_pencil_formulas(Vs)
    Vs = round(Vs, 3);
    pencil_eqs = cell(1,5);
    for i=1:5
        Ue = Vs(1:4,i);
        Ve = Vs(5:8,i);
        U_mat = sprintf('\\begin{bmatrix} %0.3f & %0.3f \\\\ %0.3f & %0.3f \\end{bmatrix}', Ue(1), Ue(3), Ue(2), Ue(4));
        V_mat = sprintf('\\begin{bmatrix} %0.3f & %0.3f \\\\ %0.3f & %0.3f \\end{bmatrix}', Ve(1), Ve(3), Ve(2), Ve(4));
        pencil_eqs{i} = sprintf('R_%d = U_%d - \\lambda V_%d &= %s - \\lambda %s', i, i, i, U_mat, V_mat);
    end
    latex_output = sprintf(['\\begin{align*}\n%s,\\\\\n%s,\\\\\n%s,\\\\\n%s,\\\\\n%s\n\\end{align*}\n'], ...
        pencil_eqs{1}, pencil_eqs{2}, pencil_eqs{3}, pencil_eqs{4}, pencil_eqs{5});
end