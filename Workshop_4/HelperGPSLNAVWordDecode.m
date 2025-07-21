function [decoded,parityBits,parityCheck] = HelperGPSLNAVWordDecode(codeword,wordNumber,lastBits)
%HelperGPSLNAVWordDecode Decode the GPS LNAV word
%
%   Note: This is a helper and its API and/or functionality may change
%   in subsequent releases.
%
%   [DECODED,PARITYBITS,PARITYCHECK] = ...
%   HelperGPSLNAVWordDecode(CODEWORD,WORDNUMBER,LASTBITS) performs error
%   detection on the CODEWORD and returns DECODED which is the first 24
%   bits of CODEWORD XORed with second element of LASTBITS. Input and
%   outputs of this function are given below:
%
%   Inputs:
%   CODEWORD - 30 bit word in subframe of GPS LNAV data. A column vector of
%              length 30.
%   WORDNUMBER - Word number (1 to 10)
%   LASTBITS - Last two bits of the previous word
%
%   Outputs:
%   DECODED - Column vector of length 24 containing decoded bits
%   PARITYCHECK - A logical value indicating if whether parity check is
%   passed. TRUE indicates parity check is passed.

%   Copyright 2021-2022 The MathWorks, Inc.

codeword = codeword(:);
decoded = double(xor(codeword(1:24),lastBits(2)));
reEncoded = HelperGPSLNavWordEnc(decoded,wordNumber,lastBits);
parityBits = reEncoded(25:end);
parityCheck = isequal(codeword(25:end),parityBits);
end

function codeWord = HelperGPSLNavWordEnc(data,wordNumber,lastBits)
%Encode the GPS LNAV word
%   CODEWORD = gpsLNavWordEnc(DATA,WORDNUMBER) encodes DATA using the
%   encoding scheme specified in IS-GPS-200K to obtain the encoded
%   CODEWORD. DATA is a vector of bits. If WORDNUMBER is either 2 or 10,
%   then DATA should be a vector of 22 bits and for other value of
%   WORDNUMBER, DATA is a vector of 24 bits. CODEWORD is a vector of 30
%   bits.

d29Star = lastBits(1);
d30Star = lastBits(2);
% Parity Polynomials for bits 25 to 30
d25 = [1 2 3 5 6 10 11 12 13 14 17 18 20 23];
d26 = [2 3 4 6 7 11 12 13 14 15 18 19 21 24];
d27 = [1 3 4 5 7 8 12 13 14 15 16 19 20 22];
d28 = [2 4 5 6 8 9 13 14 15 16 17 20 21 23];
d29 = [1 3 5 6 7 9 10 14 15 16 17 18 21 22 24];
d30 = [3 5 6 8 9 10 11 13 15 19 22 23 24];

data = data(:); % As data should always be a column vector

if wordNumber~=2 && wordNumber~=10
    %Vectorize parity application to data
    dataParity = d30Star*ones(24, 1);

    %Generate Codeword based on the parity equations from the ICD -
    %Hamming(32,26)
    codeWord = mod([dataParity+data; d29Star+sum(data(d25));...
        d30Star+sum(data(d26)); d29Star+sum(data(d27));...
        d30Star+sum(data(d28)); d30Star+sum(data(d29));...
        d29Star+sum(data(d30))],2);

end

%Parity for the Handover word and word number 10 is generated differently
% Word number 2 and 10 have bits 29 and 30 set to 0 and the parity is then
% generated accordingly. The values of bits 23 and 24 are then calculated
% to satisfy this condition.
if wordNumber==2 || wordNumber==10
    %Vectorize Parity application to data
    dataParity = d30Star*ones(22, 1);

    %Bits 23 and 24 have to be calculated by setting D29 and D30 to zero
    bit24 = mod(sum(data(d29(d29<23)))+...
        d30Star,2);
    bit23 = mod(sum(data(d30(d30<23)))+bit24+...
        d29Star,2);

    %Generate Parity Encoded Data
    codeWord = mod([dataParity+data(1:22); d30Star+bit23;...
        d30Star+bit24; d29Star+sum(data(d25(d25<23)))+bit23;...
        d30Star+sum(data(d26(d26<23)))+bit24; d29Star+sum(data(d27(d27<23)));...
        d30Star+sum(data(d28(d28<23)))+bit23; 0; 0],2);

end
end