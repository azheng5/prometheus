function sol =day2sol(y,m,d)

sol_0_date = datetime(2029, 3, 3);  % Sol 0 -> March 3, 2029
sol_670_date = datetime(2031, 1, 19); % Sol 670 -> January 19, 2031

% Calculate days per sol
days_per_sol = days(sol_670_date - sol_0_date) / 670;

% Input date (Change this to any desired date)
input_date = datetime(y, m, d); % Example: August 27, 2030

% Convert date to sol number
sol = (days(input_date - sol_0_date)) / days_per_sol;

end