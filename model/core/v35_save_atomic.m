function v35_save_atomic(f, S)
%V35_SAVE_ATOMIC Save the fields of struct S to file f atomically
% (write to a unique temporary file in the same folder, then rename).
[d, ~] = fileparts(f);
if isempty(d), d = pwd; end
if ~exist(d, 'dir'), mkdir(d); end
tmp = [tempname(d) '.mat'];
save(tmp, '-v7', '-struct', 'S');
movefile(tmp, f, 'f');
end
