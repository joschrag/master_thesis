function [hashStr] = str_to_hash(text)
%STR_TO_HASH Create a hash from given text.
if ~isa(text,"string")
    text = string(text);
end
sha256hasher = System.Security.Cryptography.SHA256Managed;
sha256 = uint8(sha256hasher.ComputeHash(uint8(char(text))));
hash_hex = dec2hex(sha256);
hashStr = string([]);
nBytes = length(hash_hex);
for k=1:nBytes
    hashStr(end+1:end+2) = hash_hex(k,:);
end
hashStr = join(hashStr,"");
end

