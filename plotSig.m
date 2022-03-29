%Getting the input M from the user

prompt = 'please enter the matrix S';
S = input(prompt)
plotSignal(S)
plotSignal(gramSchmidt(S))


%plot function
function graphSig = plotSignal(S)
	[m ,n] = size(S);

	%dividing the period to three time intervals
	x = zeros(size(m, 1));
	t1 = 0:0.1:1;
	t2 = 1:0.1:2;
	t3 = 2:0.1:3;
	t = [t1 t2 t3];
	x1 = zeros(size(t1));
	x2 = zeros(size(t2));
	x3 = zeros(size(t3));
	signal = zeros(size(t));
	

	%checking the column index
	cindex = 1;
	while cindex < n+1
		y = num2cell(x);
		%checking the row index
		for rindex = 1 : m + 1
			if rindex < m + 1
				y{rindex} = S(rindex,cindex);
			else
				x = cell2mat(y);
				x1 = x(1,1) * ones(size(t1));
				x2 = x(1,2) * ones(size(t2));
				x3 = x(1,3) * ones(size(t3));
				signal = [x1 x2 x3];
			end
		end
		figure
		plot(t, signal)
		hold on
		cindex = cindex + 1;
	end
end

% S is input matrix of M signals
% g is a matrix of n column vectors that form an orthogonal basis 
% V is output matrix (phi) of n column vectors that form an orthonormal basis
function V = gramSchmidt (S) 
	[m,n]=size(S); % m is number of rows, n is number of columns 
	g = zeros (m, n);  
	V = zeros (m, n);
  
	g (:, 1) = S (:, 1); % g1 = S1
	V (:, 1) = g (:, 1) ./ norm (g (:, 1));  % phi 1 (t) = g1 / |g1| = s1 / sqrt(E1)

	for i = 2 : n
  		sum = zeros (m, 1); % Sum of vector projections 
  		for j = 1 : i-1
    		sum = sum + coProj (S (:, i), g (:, j));
  		end
  		g (:, i) = S (:, i) - sum; % gi(t) = Si - sum of vector projections 
  		V (:, i) = g (:, i) / norm (g (:, i)); % phi i (t) = gi / |gi| = si / sqrt(Ei)
	end
end

% Define a function for obtaining the coefficients of Si(t) with orthonormal
% projections 
function C = coProj (S, g) 
[m,n]=size(S); % m is number of rows, n is number of columns 
rec = zeros (m, n); 

for i = 1 : n
  rec(:, i) = dot (g(:, i), S) ./ dot(g(:, i), g(:, i)) .* g(:, i); % Sij * phi j
end
C = sum(rec,2); %total sum of projections 
end