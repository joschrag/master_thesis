classdef FF
    %FF Summary of this class goes here
    %   Detailed explanation goes here

    properties
        value
        order double
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
                error("Order must be prime.")
            end
            obj.order = order;
            obj.value = mod(value,order);
        end
        %mpower

        function r = plus(o1,o2)
            arguments
                o1 FF
                o2 FF
            end
            if o1.order ~= o2.order
                error("Order of summands must be the same.\n" + ...
                    "Orders are %i and %i.",o1.order,o2.order)
            end
            r = FF([o1.value] + [o2.value],o1.order);
        end
        function r = minus(o1,o2)
            arguments
                o1 FF
                o2 FF
            end
            if o1.order ~= o2.order
                error("Order of summands must be the same.\n" + ...
                    "Orders are %i and %i.",o1.order,o2.order)
            end
            r = FF([o1.value] - [o2.value],o1.order);
        end
        function r = uminus(obj)
            arguments
                obj FF
            end
            r = FF(-[obj.value],obj.order);
        end
        function r = uplus(obj)
            arguments
                obj FF
            end
            r = obj;
        end
        function r = times(o1,o2)
            if isa(o1,"FF") && isa(o2,"FF")
                if o1.order ~= o2.order
                    error("Order of summands must be the same.\n" + ...
                        "Orders are %i and %i.",o1.order,o2.order)
                end
                dot_mul = [o1.value].*[o2.value];
                ord = o1.order;
            elseif isa(o1,"FF")
                dot_mul = [o1.value].*o2;
                ord = o1.order;
            else
                dot_mul = o1.*[o2.value];
                ord = o2.order;
            end
            r = FF(dot_mul,ord);
        end
        function r = mtimes(o1,o2)
            arguments
                o1 FF
                o2 FF
            end
            if o1.order ~= o2.order
                error("Order of summands must be the same.\n" + ...
                    "Orders are %i and %i.",o1.order,o2.order)
            end
            r = FF(o1.value*o2.value,o1.order);
        end
        function r = rdivide(o1,o2)
            arguments
                o1 FF
                o2 FF
            end
            if o1.order ~= o2.order
                error("Order of summands must be the same.\n" + ...
                    "Orders are %i and %i.",o1.order,o2.order)
            end
            if isa(o1,"FF") && isa(o2,"FF")
                if o1.order ~= o2.order
                    error("Order of summands must be the same.\n" + ...
                        "Orders are %i and %i.",o1.order,o2.order)
                end
                dot_div = [o1.value].\[o2.value];
                ord = o1.order;
            elseif isa(o1,"FF")
                dot_div = [o1.value].\o2;
                ord = o1.order;
            else
                dot_div = o1.\[o2.value];
                ord = o2.order;
            end
            r = FF(dot_div,ord);

        end
        function r = power(obj,pow)
            arguments
                obj FF
                pow {mustBeInteger}
            end
            if pow >= 0
                val = obj.value.^pow;
            elseif (obj.value == 1) || (obj.value == (obj.order-1))
                val = a;
            elseif ~all(obj.value)
                error('FF:power:inverse','zero has no inverse')
            else
                [~, c, ~] = gcd(1:obj.order-1, obj.order);
                y = mod(c, obj.order)';

                val = y(obj.value).^(-pow);
            end
            r = FF(val,obj.order);
        end
        function r = mldivide(o1,o2)
            arguments
                o1
                o2
            end
            if ~isa(o1,"FF")
                o1 = FF(o1,o2.order);
            end
            if ~isa(o2,"FF")
                o2 = FF(o2,o1.order);
            end
            if o1.order ~= o2.order
                error("Order of summands must be the same.\n" + ...
                    "Orders are %i and %i.",o1.order,o2.order)
            end
            %Code adapted from Matlab Comm Toolbox gflineq function
            coder.internal.errorIf(size(o1.value, 1) ~= size(o2.value, 1), 'FF:mldivide:InvalidAB');
            aa = [o1.value o2.value];
            [m_aa, n_aa] = size(aa);

            row_idx = 1;
            column_idx = 1;
            row_store = [];
            column_store = [];

            % Find the multiplicative inverse of the field elements.
            % This will be used for setting major elements in the matrix to one.
            [~, c, ~] = gcd(1:o1.order-1, o1.order);
            field_inv = mod(c, o1.order)';

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
                        aa(row_idx,:) = mod( field_inv( aa(row_idx,column_idx) ) * aa(row_idx,:), o1.order );
                    end

                    % Find the other elements in the column that must be cleared,
                    idx = find(aa(:,column_idx)~=0)';
                    % and set those elements to zero.
                    for i = idx
                        if i ~= row_idx
                            aa(i,:) = mod( aa(i,:) + aa(row_idx,:) * (o1.order - aa(i,column_idx)), o1.order );
                        end
                    end

                    column_idx = column_idx + 1;

                end

                row_idx = row_idx + 1;

            end

            if ( rank(aa) > rank(aa(:,1:size(o1.value, 2))))
                % The case of no solution.
                warning(message('FF:mldivide:NoSolution'));
                x = [];
            else
                x = zeros(size(o1.value, 2), 1);
                x(column_store,1) = aa(row_store,n_aa);
            end
            r = FF(x,o1.order);
        end
        function r = mrdivide(o1,o2)
            arguments
                o1
                o2
            end
            if ~isa(o1,"FF")
                o1 = FF(o1,o2.order);
            end
            if ~isa(o2,"FF")
                o2 = FF(o2,o1.order);
            end
            if o1.order ~= o2.order
                error("Order of summands must be the same.\n" + ...
                    "Orders are %i and %i.",o1.order,o2.order)
            end
            r = mldivide(o2',o1')';
        end
        function obj = ctranspose(obj) % overload conjugate transpose: '
            % ctranspose and transpose have identical behavior
            obj.value=obj.value';
        end
        function obj = transpose(obj)            
            obj.value=obj.value.';
        end
        function r = eq(o1,o2)
            arguments
                o1 FF
                o2 FF
            end
            r = false;
            if (o1.order == o2.order) && all(o1.value == o2.value)
                r = true;
            end
        end
        function r = ne(o1,o2)
            arguments
                o1 FF
                o2 FF
            end
            r = ~eq(o1,o2);
        end
        function r = horzcat(o1,o2)
            arguments
                o1 FF
            end
            arguments (Repeating)
                o2 FF
            end
            for o = o2
                if o{1}.order == o1.order
                    o1.value = [o1.value, o{1}.value];
                else
                    error("FF:horzcat:diffOrders")
                end
            end
            r = o1;
        end
        function r = vertcat(o1,o2)
            arguments
                o1 FF
            end
            arguments (Repeating)
                o2 FF
            end
            for o = o2
                if o{1}.order == o1.order
                    o1.value = [o1.value; o{1}.value];
                else
                    error("FF:horzcat:diffOrders")
                end
                r = o1;
            end
        end
    end
end