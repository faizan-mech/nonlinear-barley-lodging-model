function cs = v35_find_case(id)
%V35_FIND_CASE Case definition by id (see V35_CASES).
C = v35_cases();
k = find(strcmp({C.id}, id));
if isempty(k), error('V35:case','Unknown case %s.', id); end
cs = C(k);
end
