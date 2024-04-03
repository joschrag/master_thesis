classdef FF
    %FF Summary of this class goes here
    %   Detailed explanation goes here

    properties
        value
        order double
    end

    methods (Static)
        function val = simplify_mod(val,order)
            arguments
                val (1,1)
                order (1,1) {mustBeInteger,mustBePositive}
            end
            try
                l(1) = lhs(val);
                l(2) = rhs(val);
            catch
                l = val;
            end
            if numel(l) == 1
                [c,t] = coeffs(l,symvar(l));
                val = sum(mod(c,order).*t);
            elseif numel(l) ==2
                [c1,t1] = coeffs(l(1),symvar(l(1)));
                [c2,t2] = coeffs(l(2),symvar(l(2)));
                val = sum(mod(c1,order).*t1) == sum(mod(c2,order).*t2);
            end
        end
    end
    methods
        function obj = FF(value,order)
            %FF Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                value (:,:)
                order (1,1) {mustBeInteger,mustBePositive}
            end
            if ~isprime(order)
                error("FF:order","Order must be prime.")
            end
            obj.order = order;
            if isa(value,"numeric")
                obj.value = mod(value,order);
            elseif isa(value,"sym")
                obj.value = arrayfun(@(x) FF.simplify_mod(x,order),value);
            end
        end
        function r = plus(self,other)
            arguments
                self FF
                other FF
            end
            if self.order ~= other.order
                error("Order of summands must be the same.\n" + ...
                    "Orders are %i and %i.",self.order,other.order)
            end
            r = FF([self.value] + [other.value],self.order);
        end
        function r = minus(self,other)
            arguments
                self FF
                other FF
            end
            if self.order ~= other.order
                error("Order of summands must be the same.\n" + ...
                    "Orders are %i and %i.",self.order,other.order)
            end
            r = FF([self.value] - [other.value],self.order);
        end
        function r = uminus(self)
            arguments
                self FF
            end
            r = FF(-[self.value],self.order);
        end
        function r = uplus(self)
            arguments
                self FF
            end
            r = self;
        end
        function r = times(self,other)
            if isa(self,"FF") && isa(other,"FF")
                if self.order ~= other.order
                    error("Order of summands must be the same.\n" + ...
                        "Orders are %i and %i.",self.order,other.order)
                end
                dot_mul = [self.value].*[other.value];
                ord = self.order;
            elseif isa(self,"FF")
                dot_mul = [self.value].*other;
                ord = self.order;
            else
                dot_mul = self.*[other.value];
                ord = other.order;
            end
            r = FF(dot_mul,ord);
        end
        function r = mtimes(self,other)
            arguments
                self FF
                other FF
            end
            if self.order ~= other.order
                error("Order of summands must be the same.\n" + ...
                    "Orders are %i and %i.",self.order,other.order)
            end
            r = FF(self.value*other.value,self.order);
        end
        function r = rdivide(self,other)
            arguments
                self FF
                other FF
            end
            if self.order ~= other.order
                error("Order of summands must be the same.\n" + ...
                    "Orders are %i and %i.",self.order,other.order)
            end
            if isa(self,"FF") && isa(other,"FF")
                if self.order ~= other.order
                    error("Order of summands must be the same.\n" + ...
                        "Orders are %i and %i.",self.order,other.order)
                end
                dot_div = [self.value].\[other.value];
                ord = self.order;
            elseif isa(self,"FF")
                dot_div = [self.value].\other;
                ord = self.order;
            else
                dot_div = self.\[other.value];
                ord = other.order;
            end
            r = FF(dot_div,ord);

        end
        function r = power(self,pow)
            arguments
                self FF
                pow {mustBeInteger}
            end
            if pow >= 0
                val = self.value.^pow;
            elseif ~all(self.value)
                error('FF:power:inverse','zero has no inverse')
            else
                [~, c, ~] = gcd(1:self.order-1, self.order);
                y = mod(c, self.order)';

                val = y(self.value).^(-pow);
            end
            r = FF(val,self.order);
        end
        function r = mldivide(self,other)
            arguments
                self
                other
            end
            if ~isa(self,"FF")
                self = FF(self,other.order);
            end
            if ~isa(other,"FF")
                other = FF(other,self.order);
            end
            if self.order ~= other.order
                error("Order of summands must be the same.\n" + ...
                    "Orders are %i and %i.",self.order,other.order)
            end
            %Code adapted from Matlab Comm Toolbox gflineq function
            coder.internal.errorIf(size(self.value, 1) ~= size(other.value, 1), 'FF:mldivide:InvalidAB');
            aa = [self.value other.value];
            [m_aa, n_aa] = size(aa);

            row_idx = 1;
            column_idx = 1;
            row_store = [];
            column_store = [];

            % Find the multiplicative inverse of the field elements.
            % This will be used for setting major elements in the matrix to one.
            [~, c, ~] = gcd(1:self.order-1, self.order);
            field_inv = mod(c, self.order)';

            % Search for major elements, trying to form 'identity' matrix.
            while (row_idx <= m_aa) && (column_idx < n_aa)

                % Look for a major element in the current column.
                while (aa(row_idx,column_idx) == 0) && (column_idx < n_aa)

                    % In the current column, search below all the rows that already
                    % have major elements.
                    idx = find( aa(row_idx:m_aa, column_idx) ~= 0 );

                    if isempty(idx)
                        % There are no major elements in this column.
                        % Move to the next column.
                        column_idx = column_idx + 1;

                    else
                        % There is a major element in this column.
                        % See if any are already equal to one.
                        idx = [ find(aa(row_idx:m_aa, column_idx) == 1); idx ];
                        idx = idx(1);

                        % Swap the current row with a row containing a major element.
                        temp_row = aa(row_idx,:);
                        aa(row_idx,:) = aa(row_idx+idx-1,:);
                        aa(row_idx+idx-1,:) = temp_row;

                    end
                end


                % Clear all non-zero elements in the column except the major element,
                % and set the major element to one.
                if ( ( aa(row_idx,column_idx) ~= 0 ) && ( column_idx < n_aa ) )

                    % The current element is a major element.
                    row_store = [row_store row_idx];
                    column_store = [column_store column_idx];

                    % If the major element is not already one, set it to one.
                    if (aa(row_idx,column_idx) ~= 1)
                        aa(row_idx,:) = mod( field_inv( aa(row_idx,column_idx) ) * aa(row_idx,:), self.order );
                    end

                    % Find the other elements in the column that must be cleared,
                    idx = find(aa(:,column_idx)~=0)';
                    % and set those elements to zero.
                    for i = idx
                        if i ~= row_idx
                            aa(i,:) = mod( aa(i,:) + aa(row_idx,:) * (self.order - aa(i,column_idx)), self.order );
                        end
                    end

                    column_idx = column_idx + 1;

                end

                row_idx = row_idx + 1;

            end
            if ( rank(aa) > rank(aa(:,1:size(self.value, 2))))
                % The case of no solution.
                warning(message('FF:mldivide:NoSolution'));
                x = [];
            elseif isa(aa,"numeric")
                x = zeros(size(self.value, 2), 1);
                x(column_store,1) = aa(row_store,n_aa);
            else
                x = sym("x",[size(self.value, 2), 1]);
                x(column_store,1) = aa(row_store,n_aa);
            end
            r = FF(x,self.order);
        end
        function r = mrdivide(self,other)
            arguments
                self
                other
            end
            if ~isa(self,"FF")
                self = FF(self,other.order);
            end
            if ~isa(other,"FF")
                other = FF(other,self.order);
            end
            if self.order ~= other.order
                error("Order of summands must be the same.\n" + ...
                    "Orders are %i and %i.",self.order,other.order)
            end
            r = mldivide(other',self')';
        end
        function self = ctranspose(self) % overload conjugate transpose: '
            % ctranspose and transpose have identical behavior
            self.value=self.value';
        end
        function self = transpose(self)
            self.value=self.value.';
        end
        function ret = eq(self,other)
            arguments
                self FF
                other FF
            end
            if self.order ~= other.order
                error("Order of summands must be the same.\n" + ...
                    "Orders are %i and %i.",self.order,other.order)
            end
            if isa(self.value,"sym")
                s_val = FF.simplify_mod(self.value,self.order);
                o_val = FF.simplify_mod(other.value,other.order);
                ret = FF(s_val == o_val,self.order);
            elseif isa(self.value,"numeric")
                ret = false;
                if (self.order == other.order) && all(self.value == other.value)
                    ret = true;
                end
            end
        end
        function r = ne(self,other)
            arguments
                self FF
                other FF
            end
            r = ~eq(self,other);
        end
        function ret = horzcat(self,other)
            arguments
                self FF
            end
            arguments (Repeating)
                other FF
            end
            ret = self;
            for o = other
                if o{1}.order == self.order
                    ret = builtin("horzcat",ret, o{1});
                else
                    error("FF:horzcat:diffOrders")
                end
            end
        end
        function ret = vertcat(self,other)
            arguments
                self FF
            end
            arguments (Repeating)
                other FF
            end
            ret = self;
            for o = other
                if o{1}.order == self.order
                    ret = builtin("vertcat",ret, o{1});
                else
                    error("FF:vertcat:diffOrders")
                end
            end
        end
        function ret = mpower(self,pow)
            arguments
                self FF
                pow (1,1) {mustBeInteger}
            end
            if size(self.value,1) ~= size(self.value,2)
                error("FF:mpower:NonSquareMatrix")
            end
            if pow == 0
                ret = FF(eye(size(self.value)),self.order);
            elseif pow > 0
                if pow == 1
                    ret = self;
                else
                    ret = FF(self.value^pow,self.order);
                end
            else
                ret = inv(self);
                if pow < -1
                    ret = ret^(-pow);
                end
            end
        end
        function ret = inv(self)
            arguments
                self FF
            end
            [m,n] = size(self.value);
            if m ~= n
                error("FF:inv:NonSquareMatrix")
            end
            d = FF(det(self),self.order);
            if d.value == 0
                error("FF:inv:SingularMatrix")
            end
            ret = d.^(-1) .* adj(self);
        end
        function ret = adj(self)
            arguments
                self FF
            end
            [m,n] = size(self.value);
            if m ~= n
                error("FF:adj:NonSquareMatrix")
            end
            A = zeros(m,n);
            for i=1:m
                for j=1:n
                    A(j,i) = cofactor(self,i,j);
                end
            end
            ret = FF(A,self.order);
        end
        function ret = cofactor(self,i,j)
            arguments
                self FF
                i (1,1) {mustBeInteger,mustBePositive}
                j (1,1) {mustBeInteger,mustBePositive}
            end
            [m,n] = size(self.value);
            if m ~= n
                error("FF:cofactor:NonSquareMatrix")
            end
            validateattributes(i,"numeric",{'>',0,'<=',m},"FF:cofactor")
            validateattributes(j,"numeric",{'>',0,'<=',n},"FF:cofactor")
            ret = mod(round((-1)^(i+j)*det(self.value([1:i-1 i+1:n],[1:j-1 j+1:n]))),self.order);
        end
        function r =  rank(self)
            %GFRANK Compute the rank of a matrix over a Galois field.
            %
            %   RK = GFRANK(A) calculates the rank of the matrix A in GF(2).
            %
            %   RK = GFRANK(A, P) calculates the rank of the matrix A in GF(P).
            %
            %   Note: This function performs computations in GF(P) where P is prime. To
            %   work in GF(2^M), use the RANK function with Galois arrays.
            %
            %   See also GFLINEQ, GF/RANK.

            %   The method used here is similar to the Gaussian elimination. The
            %   algorithm has taken advantage of the binary computation, and double
            %   sided elimination has been used.

            %   Copyright 1996-2017 The MathWorks, Inc.
            [n, m] = size(self.value);
            if n < m
                self.value = self.value';
            end
            p = self.order;
            x = self.value(:);
            if ((max(x) >=p) || (min(x) < 0) || any(any(floor(x)~=x)))
                error(message('comm:gfrank:InvalidElementsInAB'))
            end
            [n, m] = size(self.value);
            k = 1;
            i = 1;
            ii = [];
            kk = [];

            % forward major element selection
            while (i <= n) && (k <= m)
                % make the diagonal element into 1, or jump one line.
                while (self.value(i,k) == 0) && (k <= m)
                    ind = find(self.value(i:n, k) ~= 0);
                    if isempty(ind) && (k == m)
                        break;
                    elseif isempty(ind)
                        k = k + 1;
                    else
                        indx = find(self.value(i:n, k) == 1);
                        if isempty(indx)
                            ind_major = ind(1);
                        else
                            ind_major = indx(1);
                        end
                        j = i + ind_major - 1;
                        tmp = self.value(i, :);
                        self.value(i,:) = self.value(j, :);
                        self.value(j, :) = tmp;
                    end
                end

                % clear all nonzero elements in the column except the major element.
                if (self.value(i,k) ~= 0)
                    % to make major element into 1
                    if (self.value(i,k) ~= 1)
                        self.value(i,:) = rem(self.value(i,k)^(p-2) * self.value(i,:), p);
                    end
                    ind = find(self.value(:,k) ~= 0)';
                    ii = [ii i];
                    kk = [kk k];
                    vec = (k:m);
                    for j = ind
                        if j > i
                            % to make sure the column will be zero except the major element.
                            self.value(j, vec) = rem(self.value(j, vec) + self.value(i, vec) * (p - self.value(j, k)), p);
                        end
                    end
                    k = k + 1;
                end
                i = i + 1;
            end

            r = find(sum(self.value,2)>0,1,'last');
            if isempty(r)
                r = 0;
            end
        end
        function ret = subsref(self,S)
            for s = S
                self = builtin("subsref",self,s);
            end
            ret = self;
        end
        function res = sum(self)
            res = FF(sum(self.value),self.order);
        end
        function ret = expand(self,varargin)
            val = arrayfun(@(x) FF.simplify_mod(x,self.order),self.value);
            if numel(varargin) == 0
                ret = FF(expand(val),self.order);
            else
                ret = FF(expand(val,varargin),self.order);
            end
        end
        function ret = simplify(self,varargin)
            val = arrayfun(@(x) FF.simplify_mod(x,self.order),self.value);
            if numel(varargin) == 0
                ret = FF(simplify(val),self.order);
            else
                ret = FF(simplify(val,varargin),self.order);
            end
        end
        function ret = subs(self,old,new)
            val = arrayfun(@(x) FF.simplify_mod(x,self.order),self.value);
            ret = FF(subs(val,old,new),self.order);
        end
        function ret = collect(self,varargin)
            val = arrayfun(@(x) FF.simplify_mod(x,self.order),self.value);
            if numel(varargin) == 0
                ret = FF(collect(val),self.order);
            else
                ret = FF(collect(val,varargin),self.order);
            end
        end
        function ret = det(self)
            if isa(self.value,"sym")
                ret = self.simplify_mod(det(self.value),self.order);
            else
                ret = mod(det(self.value),self.order);
                if (round(ret)~=ret)
                    if (abs(round(ret)-ret) < 10^-8)
                        ret = round(ret);
                    else
                        error("FF:det","Result is not an integer.")
                    end
                end
            end
        end
    end
end