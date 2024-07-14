classdef FF
    %FF Summary of this class goes here
    %   Detailed explanation goes here

    properties
        value
        order double
        vars
    end

    methods (Static)
        function val = simplify_mod(val,order,vars)
            arguments
                val (1,1)
                order (1,1) {mustBeInteger,mustBePositive}
                vars (1,:)
            end
            [c,t] = coeffs(val,vars);
            val = sum(mod(c,order).*t);
        end
        function is_same = check_arr_order(arr)
            arguments
                arr
            end
            is_same = all(arr(:)==arr(1));
        end
        function check_args(self,other)
            self_order = self.get_order();
            other_order = other.get_order();
            if ~FF.check_arr_order(self_order)
                error("FF:check_args","Order of 'self' array elements differ!")
            end
            if ~FF.check_arr_order(other_order)
                error("FF:check_args","Order of 'other' array elements differ!")
            end

            if unique(unique(self_order)) ~= unique(unique(other_order))
                error("FF:check_args","Order of arguments must be the same.\n" + ...
                    "Orders are %i and %i.",self.order,other.order)
            end
        end
    end
    methods
        function obj = FF(value,order,vars)
            %FF Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                value (:,:)
                order (1,1) {mustBeInteger,mustBePositive}
                vars (1,:) = [];
            end
            if ~isprime(order)
                error("FF:order","Order must be prime.")
            end
            obj.order = order;
            if isa(value,"numeric")
                obj.value = mod(sym(value),order);
            elseif isa(value,"sym")
                obj.vars = symvar(value);
                if numel(obj.vars) > 0
                    obj.value = arrayfun(@(x) FF.simplify_mod(x,order,obj.vars),value);
                else
                    obj.value = mod(value,order);
                    obj.vars = [];
                end
            end
        end

        function ret = get_value(self)
            arguments
                self FF
            end
            ret = arrayfun(@(x) x.value,self,UniformOutput=false);
            ret = [ret{:}];
        end

        function ret = get_order(self)
            arguments
                self FF
            end
            ret = arrayfun(@(x) x.order,self);
        end

        function ret = to_matrix(self)
            arguments
                self FF
            end
            ret = arrayfun(@(x) FF(x,self.order),self.value);
        end

        function ret = from_matrix(self)
            arguments
                self FF
            end
            arr_order = unique(arrayfun(@(x) x.order,self));
            assert(numel(arr_order)==1,"Different orders in the array!")
            ret = FF(arrayfun(@(x) x.value,self),arr_order);
        end

        function r = plus(self,other)
            arguments
                self FF
                other FF
            end
            FF.check_args(self,other)
            if numel(self.vars) > 0 || numel(other.vars) > 0
                val = arrayfun(@(x) FF.simplify_mod(x,self.order,unique([self.vars,other.vars])),self.value+other.value);
            else
                val = mod(self.value+other.value,self.order);
            end
            r = FF(val,self.order);
        end
        function r = minus(self,other)
            arguments
                self FF
                other FF
            end
            FF.check_args(self,other)
            if numel(self.vars) > 0 || numel(other.vars) > 0
                val = arrayfun(@(x) FF.simplify_mod(x,self.order,unique([self.vars,other.vars])),self.value-other.value);
            else
                val = mod(self.value-other.value,self.order);
            end
            r = FF(val,self.order);
        end
        function r = uminus(self)
            arguments
                self FF
            end
            r = arrayfun(@(x) FF(-x.value,x.order),self);
        end
        function r = uplus(self)
            arguments
                self FF
            end
            r = self;
        end
        function r = times(self,other)
            if isa(self,"FF") && isa(other,"FF")
                FF.check_args(self,other)
                dot_mul = [self.value].*[other.value];
                ord = self.order;
            elseif isa(self,"FF")
                dot_mul = [self.value].*other;
                ord = self.order;
            else
                dot_mul = self.*[other.value];
                ord = other.order;
            end
            r = FF(dot_mul,ord,unique([self.vars,other.vars]));
        end
        function self = mtimes(self,other)
            arguments
                self FF
                other FF
            end
            if self.order ~= other.order
                error("Order of summands must be the same.\n" + ...
                    "Orders are %i and %i.",self.order,other.order)
            end
            self.value = self.value*other.value;
            if isa(self.value,"sym")
                if numel([self.vars,other.vars]) > 0
                    self.value = arrayfun(@(x) FF.simplify_mod(x,self.order,unique([self.vars,other.vars])),self.value);
                else
                    self.value = mod(self.value,self.order);
                end
            else
                self.value = mod(self.value,self.order);
            end
            self.vars = unique([self.vars,other.vars]);
        end
        function r = rdivide(self,other)
            arguments
                self
                other
            end
            if isa(self,"FF") && isa(other,"FF")
                if self.order ~= other.order
                    error("Order of summands must be the same.\n" + ...
                        "Orders are %i and %i.",self.order,other.order)
                end
                other_inv = other.^(-1);
                dot_div = [self.value].*[other_inv.value];
                ord = self.order;
                n_vars = unique([self.vars,other.vars]);
            elseif isa(self,"FF")
                other_inv = FF(other.^(-1),self.order);
                dot_div = [self.value].*[other_inv.value];
                ord = self.order;
                n_vars = self.vars;
            else
                other_inv = other.^(-1);
                dot_div = self.*[other_inv.value];
                ord = other.order;
                n_vars = other.vars;
            end
            r = FF(dot_div,ord,n_vars);

        end
        function r = ldivide(self,other)
            arguments
                self
                other
            end
            r = other ./ self;

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
            r = FF(val,self.order,self.vars);
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
                        if ~isa(aa,"sym")
                            aa(row_idx,:) = mod( field_inv( aa(row_idx,column_idx) ) * aa(row_idx,:), self.order );
                        else
                            if numel(self.vars) > 0
                                aa(row_idx,:) = arrayfun(@(x) FF.simplify_mod(x,self.order,unique([self.vars,other.vars])),field_inv( aa(row_idx,column_idx) ) * aa(row_idx,:));
                            else
                                aa(row_idx,:) = mod( field_inv( aa(row_idx,column_idx) ) * aa(row_idx,:), self.order );
                            end
                        end
                    end

                    % Find the other elements in the column that must be cleared,
                    idx = find(aa(:,column_idx)~=0)';
                    % and set those elements to zero.
                    for i = idx
                        if i ~= row_idx
                            if ~isa(aa,"sym")
                                aa(i,:) = mod( aa(i,:) + aa(row_idx,:) * (self.order - aa(i,column_idx)), self.order );

                            else
                                if numel(self.vars) > 0
                                    aa(i,:) = arrayfun(@(x) FF.simplify_mod(x,self.order,unique([self.vars,other.vars])),aa(i,:) + aa(row_idx,:) * (self.order - aa(i,column_idx)));
                                else
                                    aa(i,:) = mod( aa(i,:) + aa(row_idx,:) * (self.order - aa(i,column_idx)), self.order );
                                end
                            end

                        end
                    end

                    column_idx = column_idx + 1;

                end

                row_idx = row_idx + 1;

            end
            if ( rank(aa) > rank(aa(:,1:size(self.value, 2))))
                % The case of no solution.
                warning("FF:mldivide","NoSolution");
                x = [];
            elseif isa(aa,"numeric")
                x = zeros(size(self.value, 1),size(other.value, 2));
                x(column_store,1:size(other.value, 2)) = aa(row_store,n_aa-size(other.value, 2)+1:n_aa);
            else
                x = sym("x",[size(self.value, 1),size(other.value, 2)]);
                x(column_store,1:size(other.value, 2)) = aa(row_store,n_aa-size(other.value, 2)+1:n_aa);
            end
            r = FF(x,self.order,unique([self.vars,other.vars]));
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
            if self.get_order() ~= other.get_order()
                ret = false;
                return
            end
            if size(self.value) == size(other.value)
                if numel(self.vars) > 0
                    s_val = arrayfun(@(x) FF.simplify_mod(x,self.order,self.vars),self.value);
                else
                    s_val = mod(self.value,self.order);
                end
                if numel(other.vars) > 0
                    o_val = arrayfun(@(x) FF.simplify_mod(x,other.order,other.vars),other.value);
                else
                    o_val = mod(other.value,other.order);
                end
                ret = all(all(isAlways(s_val==o_val)));

            else
                ret = false;
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
                error("FF:mpower:NonSquareMatrix","Matrix must be square.")
            end
            if pow == 0
                ret = FF(eye(size(self.value)),self.order,self.vars);
            elseif pow > 0
                if pow == 1
                    ret = self;
                else
                    ret = FF(self.value^pow,self.order,self.vars);
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
                error("FF:inv:NonSquareMatrix","Matrix must be square.")
            end
            d = FF(det(self),self.order);
            if d.value == 0
                error("FF:inv:SingularMatrix","Matrix is singular.")
            end
            ret = d.^(-1) .* adj(self);
        end
        function ret = adj(self)
            arguments
                self FF
            end
            [m,n] = size(self.value);
            if m ~= n
                error("FF:adj:NonSquareMatrix","Matrix must be square.")
            end
            A = zeros(m,n);
            for i=1:m
                for j=1:n
                    A(j,i) = cofactor(self,i,j);
                end
            end
            ret = FF(A,self.order,self.vars);
        end
        function ret = cofactor(self,i,j)
            arguments
                self FF
                i (1,1) {mustBeInteger,mustBePositive}
                j (1,1) {mustBeInteger,mustBePositive}
            end
            [m,n] = size(self.value);
            if m ~= n
                error("FF:cofactor:NonSquareMatrix","Matrix must be square.")
            end
            validateattributes(i,"numeric",{'>',0,'<=',m},"FF:cofactor")
            validateattributes(j,"numeric",{'>',0,'<=',n},"FF:cofactor")
            d = det(self.value([1:i-1 i+1:n],[1:j-1 j+1:n]));
            if (round(d)~=d)
                if (abs(round(d)-d) < 10^-8)
                    d = round(d);
                else
                    error("FF:cofactor","Subdeterminant is not an integer.")
                end
            end
            ret = mod((-1)^(i+j)*d,self.order);
        end
        function r =  rank(self)
            [n, m] = size(self.value);
            if n < m
                self.value = self.value';
            end
            p = self.order;
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
            if numel(self.vars) > 0
                val = arrayfun(@(x) FF.simplify_mod(x,self.order,self.vars),self.value);
            else
                val = mod(self.value,self.order);
            end
            if numel(varargin) == 0
                ret = FF(expand(val),self.order);
            else
                ret = FF(expand(val,varargin),self.order,self.vars);
            end
        end
        function ret = simplify(self,varargin)
            if numel(self.vars) > 0
                val = arrayfun(@(x) FF.simplify_mod(x,self.order,self.vars),self.value);
            else
                val = mod(self.value,self.order);
            end
            if numel(varargin) == 0
                ret = FF(simplify(val),self.order);
            else
                ret = FF(simplify(val,varargin),self.order,self.vars);
            end
        end
        function self = subs(self,old,new)
            if isa(new,"sym")
                self.value = FF.simplify_mod(subs(self.value,old,new),self.order,unique([self.vars,symvar(new)]));
            else
                self.value = mod(subs(self.value,old,new),self.order);
            end
        end
        function ret = collect(self,varargin)
            val = arrayfun(@(x) FF.simplify_mod(x,self.order,self.vars),self.value);
            if numel(varargin) == 0
                ret = FF(collect(val),self.order,self.vars);
            else
                ret = FF(collect(val,varargin),self.order,self.vars);
            end
        end
        function ret = det(self)
            if isa(self.value,"sym")
                ret = self.simplify_mod(det(self.value),self.order,self.vars);
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
        function ret = abs(self)
            arguments
                self FF
            end
            ret = arrayfun(@(x) FF(abs(x.value),x.order),self);
        end
        function A = rref(self)
            i = 1;
            j = 1;
            jb = zeros(1,0);
            if any(size(self)>1)
                self = self.from_matrix();
                A = self.value;
            elseif any(size(self.value)>1)
                A = self.value;
            end
            [m,n] = size(A);
            while i <= m && j <= n
                % Find value and index of largest element in the remainder of column j.
                [p, k] = max(abs(A(i:m,j)));
                k = k+i-1;
                if p == 0
                    % The column is negligible, zero it out.
                    A(i:m,j) = 0;
                    j = j + 1;
                else
                    % Remember column index
                    jb = [jb j]; %#ok<AGROW>
                    % Swap i-th and k-th rows.
                    A([i k],j:n) = A([k i],j:n);
                    % Divide the pivot row by the pivot element.
                    tmp = FF(A(i,j:n),self.order)./FF(A(i,j),self.order);
                    A(i,j:n) = tmp.value;
                    % Subtract multiples of the pivot row from all the other rows.
                    for k = [1:i-1 i+1:m]
                        A(k,j:n) = mod(A(k,j:n) - A(k,j).*A(i,j:n),self.order);
                    end
                    i = i + 1;
                    j = j + 1;
                end
            end
        end
    end
end