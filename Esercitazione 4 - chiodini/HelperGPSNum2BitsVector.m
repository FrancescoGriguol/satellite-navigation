function bits = HelperGPSNum2BitsVector(val, numBits, scaleFactor)
%Wrapper to scale the value and convert the value to bits using left MSB.
%
%   Note: This is a helper and its API and/or functionality may change
%   in subsequent releases.
%
%   BITS = HelperGPSNum2BitsVector(VAL, NUMBITS, SCALEFACTOR) scales VAL by
%   SCALEFACTOR and converts the scaled value which is in decimal format to
%   binary format keeping most significant bit (MSB) to the left. Number of
%   bits in the binary format is given by NUMBITS. If VAL is a negative
%   number, then the result is in 2's complement form.

%   Copyright 2020-2021 The MathWorks, Inc.

valScaled = round(val/scaleFactor);
if valScaled<0
    valScaled = 2^numBits+valScaled; % For converting into 2's compliment form
end

bits = rem(floor(valScaled*pow2(1-numBits:0)'),2);
end