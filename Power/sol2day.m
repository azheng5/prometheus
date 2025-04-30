function date =sol2day(sol)

%Sol to Date Converter
sol_0_date = datetime(2029, 3, 3);  % March 3, 2029
sol_670_date = datetime(2031, 1, 19); % January 19, 2031

% Calculate days per sol
days_per_sol = days(sol_670_date - sol_0_date) / 670;

% Compute the corresponding date
date = sol_0_date + days(sol * days_per_sol);

end

