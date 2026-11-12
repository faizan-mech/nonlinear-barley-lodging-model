function P = v35_apply_variant(P, ov)
%V35_APPLY_VARIANT Apply a case's parameter overrides to P.
% ov = {'path', value, ...}. 'a.b' sets P.a.b = value; 'a.b*' multiplies it.
% Unknown fields are an error, so a typo cannot silently leave the
% reference value in place.
if isempty(ov), return; end
if mod(numel(ov), 2) ~= 0, error('V35:variant','Overrides must be name/value pairs.'); end
for i = 1:2:numel(ov)
    path = ov{i}; val = ov{i+1};
    mult = path(end) == '*';
    if mult, path = path(1:end-1); end
    parts = strsplit(path, '.');
    cur = getpath(P, parts, path);
    if mult
        if ~isnumeric(cur) || isempty(cur)
            error('V35:variant','Cannot scale %s.', path);
        end
        val = cur*val;
    end
    P = setpath(P, parts, val);
end
end

function v = getpath(S, parts, full)
v = S;
for k = 1:numel(parts)
    if ~isstruct(v) || ~isfield(v, parts{k})
        error('V35:variant','Unknown parameter %s.', full);
    end
    v = v.(parts{k});
end
end

function S = setpath(S, parts, val)
if numel(parts) == 1
    S.(parts{1}) = val;
else
    S.(parts{1}) = setpath(S.(parts{1}), parts(2:end), val);
end
end
