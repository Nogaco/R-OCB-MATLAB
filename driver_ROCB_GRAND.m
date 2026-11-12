%% OCB_driver.m
% Test driver for OCB encryption/decryption
% Parameters: k=112, n=128, Nonce length = 128
clc; clear;

% Add GRAND folder to path
addpath('AES_update/');
% addpath('OCB/');
addpath('GRAND-MATLAB-main/');
% addpath('GRAND-MATLAB-main/CODES/');
addpath('Modes_GRAND');
addpath('ROCB-MATLAB-code/');
% addpath to save the results.
addpath('RESULTS/');

% Open file for appending
fileID = fopen('Report.txt', 'a');  % 'a' = append, 'w' = will overwrite if it exists
% Check if file opened successfully
if fileID == -1
    error('Failed to open file.');
end

% Chose the decoder from: 
% 'GRAND' (hard detection); 
% 'ORBGRAND' (soft detection); 
% 'ORBGRAND1' (soft detection);
DECODER='GRAND';

% Sim range in Eb/N0
ebn0=1:0.5:8;

% Run simulation until this many errors observed for each Eb/N0 value
err_thresh = 80;

% Modulation schemes available using MATLAB's toolbox, which will be used
% in a complex-valued channel
modlist = {'pi/2-BPSK','BPSK','QPSK','16QAM','64QAM','256QAM', 'BSC'}; 
bpsList = [1 1 2 4 6 8 1];
% Pick the modulation
modulation = 'BSC';
% Determine the number of bits per symbol in the modulation
nmodbits = bpsList(strcmpi(modlist,modulation));


% Pick the code
% ROCB, ROCB_CRC, RCBC or RCTR, RXTS. 
code_class = 'ROCB';

% Chose the padding from: ZERO, CRC
if isequal(code_class,'ROCB_CRC')
    PADDING='CRC';
else 
    PADDING='ZERO';
end

% Record the code_class in the data construct
code.class = code_class;

%GRAND parameters
max_query=inf; % Can reduce to an abandonment value.

%% Parameters
n = 128;
k = 116; % 112 or 116 or 120

% Reminder that this implementation is not optimized or parallelized, so 
% larger n-k may take time to run.
if n-k>16
    disp('GRAND and ORBGRAND are readily highly parallelizable in software and hardware.')
    disp('The MATLAB implementations here are for instruction and are NOT parallelized.')
    disp(['If GRAND algorithms find an error, they do so after ' ...
        'approximately 2^{n-k} code-book queries'])
    disp('and so 2^{n-k} is an upper bound on complexity in terms of query numbers.')
    disp(['With n-k=' num2str(n-k) ' expect slow decoding in noisy channels with this inefficient implementation.'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert Eb/N0 to SNR
snr_db = ebn0+10*log10(k/n)+10*log10(nmodbits);

% Pad modulation with meaningless bits if necessary.
mod_pad=zeros(1,mod(n,nmodbits));

%For each SNR record the following
num_decoded = zeros(1,length(snr_db)); % Number of packets decoded
num_demod_errs = num_decoded; % Number of demodulations that are in error
num_demod_bit_errs = num_decoded; % Number of erroneous demodulated bits
num_errs = num_decoded; % Number of erroneous decodings
num_bit_errs = num_decoded; % Number of erroneous bits
num_aband = num_decoded; % Number of abandoned decodings
num_queries = num_decoded; % Total number of code-book queries

c_num_bit_errs = num_decoded; % Number of erroneous decoded ciphers bits

% Key (random 128-bit in hex)
key = upper(reshape(dec2hex(uint8(randi([0,255],1,16)), 2).', 1, []));
% key = '218D3BCC20C7DA14238C3CB12452513C';

% Nonce (random 128-bit in hex)
N = upper(reshape(dec2hex(uint8(randi([0,255],1,16)), 2).', 1, []));

if isequal(PADDING, 'ZERO')
    % Padding - vector of n-k zeros
    P = zeros(1, n-k); % ZERO
elseif isequal(PADDING, 'CRC')
    %% CRC redundancy for R-OCB

    switch k
        case 112
            hex_poly = '0x9eb2';   % 16 CRC bits
        case 116
            hex_poly = '0x8f3';    % 12 CRC bits
        case 120
            hex_poly = '0x83';     % 8 CRC bits
        otherwise
            error('No CRC polynomial configured for k = %d.', k);
    end

    poly = koopman2matlab(hex_poly);

    [G_crc, H_crc, n_crc] = make_CRC_GH(k, poly);

    if n_crc ~= n
        error('CRC code length is %d, expected n = %d.', n_crc, n);
    end

    CRC.G = G_crc;
    CRC.H = H_crc;
    CRC.poly = poly;
    CRC.hex_poly = hex_poly;
    CRC.k = k;
    CRC.n = n;
else
    error('No CRC polynomial configured for k = %d.', k);
end

if isequal(code_class,'ROCB') || isequal(code_class,'ROCB_CRC')
    num_messages           = zeros(1,length(snr_db));
    num_rejected           = zeros(1,length(snr_db));
    num_correct_accepted   = zeros(1,length(snr_db));
    num_incorrect_accepted = zeros(1,length(snr_db));
end


for ii=1:length(snr_db)

    if isequal(modulation,'BSC')
        % Noise probability
        EbN0_coded_dB = snr_db(ii) - 10*log10(nmodbits);
        p = qfunc( sqrt( 2*10.^(EbN0_coded_dB/10) ) );
    else
        % Noise variance
        sigma2 = 1/(10^(0.1*snr_db(ii))); 
        % Using MATLAB's channel function
        awgnchan = comm.AWGNChannel('NoiseMethod','Variance','Variance',sigma2);
    end

    % Keep decoding pacekts until you see err_thresh erroneous decodings
    while num_errs(ii)<err_thresh
        % Parameters: n=128, k=112/116/120, Nonce length = 128
        % Generate random message of length between (k) and (40k) bits
        min_len = 0;
        max_len = 20*k;
        mlen_bits = randi([min_len, max_len]);
%         mlen_bits = randi([1, 30]) * k + randi([0, k-1]);
        % mlen_bits = 5*k+(k-4);
        
        % Create random binary message
        M = randi([0, 1], 1, mlen_bits, 'uint8');
        
        % Tag length (n bits, can be adjusted)
        taglen = n;
        
        % fprintf('Testing OCB with message length: %d bits\n', mlen_bits);
        
        %% Encryption
        
        if isequal(code_class,'ROCB')
            [C_out, T_out, info] = rocb_encrypt_v6(key, N, M, taglen, n, k, P);
        elseif isequal(code_class,'ROCB_CRC')
            [C_out, T_out, info] = rocb_encrypt_CRC_v1(key,N,M,taglen,n,k,CRC);
        elseif isequal(code_class,'RCBC')
            [C_out] = rcbc_encrypt_v1(key, N, M, n, k);
        elseif isequal(code_class,'RCTR')
            [C_out] = rctr_encrypt_v1(key, M, n, k);
        elseif isequal(code_class,'RXTS')
            % Key (random 128-bit in hex)
            K2_hex = upper(reshape(dec2hex(uint8(randi([0,255],1,16)), 2).', 1, []));
            sectorNum = 1;
            [C_out] = xts_encrypt_v2(key, K2_hex, M, sectorNum, n, k);
        end
        
        C =  C_out;
        T = T_out;

        % R-OCB uses taglen = n in the current simulations, so the tag is one
        % complete n-bit block. Append it to the ciphertext for channel simulation.
        TX = [C T];

        % Number of ciphertext/data blocks (tag excluded)
        m_blocks = ceil(length(C_out)/n);
        % Total number of transmitted n-bit blocks, including the tag
        tx_blocks = ceil(length(TX)/n);

        % Channel outputs
        % GRAND uses hard binary decisions.
        TX_GRAND = zeros(size(TX), 'uint8');

        % ORBGRAND uses soft information.
        TX_ORBGRAND = zeros(size(TX), 'double');
    
         for block = 1:tx_blocks
        
            idx_start = n*(block-1) + 1;
            idx_end   = min(n*block, length(TX));
        
            x = TX(idx_start:idx_end);
        
            if isequal(modulation,'BSC')
                %% BSC
                % --- Channel: flip bits with probability p ---
                % Keep c as row vector and ensure it's binary (0/1, not logical)
                x_binary = double(x > 0);  % Ensure row vector of 0s and 1s
                
                % Generate flip mask with same dimensions
                flipMask = rand(size(x_binary)) < p;
                
                % Apply channel noise (XOR with flip mask)
                rBits = xor(x_binary, flipMask);
                
                % --- Soft information for ORBGRAND (BSC LLRs) ---
                p_eff = min(max(p,1e-12), 1-1e-12);      % avoid log(0)
                LLR0  = log((1 - p_eff) / p_eff);         % constant magnitude
                y_soft = (1 - 2*rBits) * LLR0;           % +LLR0 for 0, -LLR0 for 1
    
                % Hard demodulate (ensure same type as c_binary)
                y_demod = rBits;  % Already 0/1 double array
                
                % Update c to match the format expected by comparison
                c = x_binary;
            else
                %% AWGN + BPSK
 
	            % Modulate - BPSK modulation: 0 -> +1, 1 -> -1
	            % modOut = nrSymbolModulate([c mod_pad]',modulation);
                modOut = nrSymbolModulate(double([x mod_pad]'), modulation);
                % Add White Gaussian noise
                rSig = awgnchan(modOut);
                % Soft demodulate
                y_soft = nrSymbolDemodulate(rSig,modulation,sigma2);
    	        y_soft = y_soft(1:n)';
                % Hard demodulate
                y_demod = (y_soft<0);
            end
            %% Channel-error statistics
            % Count the times the demodulation is in error
            x_logical = logical(x);

            if ~isequal(logical(y_demod), x_logical)
                num_demod_errs(ii) = num_demod_errs(ii)+1;
                % Count the demodulated bit errors
                num_demod_bit_errs(ii) = num_demod_bit_errs(ii)+sum(logical(y_demod) ~= x_logical);
            end
        
            TX_GRAND(idx_start:idx_end) = uint8(y_demod);

            TX_ORBGRAND(idx_start:idx_end) = double(y_soft);
         end

        %% Separate noisy ciphertext and noisy tag again
        
        C_len = length(C);
        T_len = length(T);
        
        % Hard-decision received values
        C_in_GRAND = TX_GRAND(1:C_len);
        T_in_GRAND = TX_GRAND(C_len+1:C_len+T_len);
        
        % Soft received values
        C_in_ORBGRAND = TX_ORBGRAND(1:C_len);
        T_in_ORBGRAND = TX_ORBGRAND(C_len+1:C_len+T_len);

        %% Decryption
        if(~isequal(length(C_in_GRAND), length(C)))
            disp('1')
        end
        % Decode with GRAND (hard detection)
        if isequal(DECODER,'GRAND')
            
            C_in = C_in_GRAND;
            if isequal(code_class,'ROCB')
                % [c_decoded,M_dec,valid,abandoned,n_guesses] = ocb_decrypt_GRAND_v4(key, N, C_in, T_in, taglen, n, k, max_query);
                [c_decoded, M_dec, valid, abandoned, n_guesses, info] = rocb_decrypt_GRAND_v6(key, N, C_in, T_in_GRAND, taglen, n, k, max_query, P);
                if ~valid
                    num_rejected(ii) = num_rejected(ii) + 1;
                
                elseif isequal(M_dec, M)
                    num_correct_accepted(ii) = num_correct_accepted(ii) + 1;
                
                else
                    num_incorrect_accepted(ii) = num_incorrect_accepted(ii) + 1;
                end
            elseif isequal(code_class,'ROCB_CRC')
                [c_decoded,M_dec,valid,abandoned,n_guesses,info] = rocb_decrypt_GRAND_CRC_v1( ...
                        key,N,C_in,T_in_GRAND,taglen,n,k,max_query,CRC);
                if ~valid
                    num_rejected(ii) = num_rejected(ii) + 1;
                elseif isequal(M_dec, M)
                    num_correct_accepted(ii) = num_correct_accepted(ii) + 1;
                else
                    num_incorrect_accepted(ii) = num_incorrect_accepted(ii) + 1;
                end
            elseif isequal(code_class,'RCBC')
                [c_decoded,M_dec,abandoned,n_guesses,info] = rcbc_decodedrypt_GRAND_v1(key, N, C_in, n, k, max_query);
            elseif isequal(code_class,'RCTR')
                [c_decoded,M_dec,abandoned,n_guesses] = rctr_decrypt_GRAND_v1(key, C_in, n, k, max_query);
            elseif isequal(code_class,'RXTS')
                [c_decoded,M_dec,abandoned,n_guesses] = xts_decrypt_GRAND_v2(key, K2_hex, C_in, sectorNum, n, k, max_query);
            end
        % Decode with basic ORBGRAND (soft detection)
        elseif isequal(DECODER,'ORBGRAND')
            C_in = C_in_ORBGRAND;
            if isequal(code_class,'ROCB')
                [c_decoded,M_dec,valid,abandoned,n_guesses,info] = rocb_decrypt_ORBGRAND_v6(key, N, C_in, T_in_ORBGRAND, taglen, n, k, max_query, P);
                if ~valid
                    num_rejected(ii) = num_rejected(ii) + 1;
                
                elseif isequal(M_dec, M)
                    num_correct_accepted(ii) = num_correct_accepted(ii) + 1;
                
                else
                    num_incorrect_accepted(ii) = num_incorrect_accepted(ii) + 1;
                end
            elseif isequal(code_class,'ROCB_CRC')
                [c_decoded,M_dec,valid,abandoned,n_guesses,info] = rocb_decrypt_ORBGRAND_CRC_v1( ...
                    key,N,C_in,T_in_ORBGRAND,taglen,n,k,max_query,CRC);
    
                if ~valid
                    num_rejected(ii) = num_rejected(ii) + 1;
                elseif isequal(M_dec, M)
                    num_correct_accepted(ii) = num_correct_accepted(ii) + 1;
                else
                    num_incorrect_accepted(ii) = num_incorrect_accepted(ii) + 1;
                end
            elseif isequal(code_class,'RCBC')
                [c_decoded,M_dec,abandoned,n_guesses] = rcbc_decodedrypt_ORBGRAND_v1(key, N, C_in, n, k, max_query);
            elseif isequal(code_class,'RCTR')
                [c_decoded,M_dec,abandoned,n_guesses] = rctr_decrypt_ORBRAND_v1(key, C_in, n, k, max_query);
            end
        else
            disp('DECODER must be GRAND or ORBGRAND or ORBGRAND1')
            return;
        end      

        if ~valid
            fprintf('ROCB invalid: %s\n', info.decode_failure_reason);
        end

        % Total number of queries made at this SNR
        num_queries(ii) = num_queries(ii) + sum(n_guesses(:));

        num_decoded(ii)=num_decoded(ii)+m_blocks;

        % Complete R-OCB message-level outcome
        num_messages(ii) = num_messages(ii) + 1;
        
        %% Validation
        if ~isequal(c_decoded, C)

            for block=1:m_blocks
                if (block < m_blocks)
                    if ~isequal(c_decoded(n*(block-1)+1:n*block),C(n*(block-1)+1:n*block))
                        num_errs(ii)=num_errs(ii)+1;
                    end
                else
                    if ~isequal(c_decoded(n*(block-1)+1:end),C(n*(block-1)+1:end))
                        num_errs(ii)=num_errs(ii)+1;
                    end
                end
            end 
            num_bit_errs(ii) = num_bit_errs(ii) + sum(double(c_decoded) ~= double(C));
            % fprintf('Decryption failed! Message integrity compromised.\n');

            % Find the positions where the bits are different
%             diff_positions = find(c_decoded ~= C); 
            % % Display the positions
            % disp('Different bit positions are:');
            % disp(diff_positions);
       
    
            % Report to the terminal every time a new error is observed
            % disp(['(end) mlen_bits = s*k : n=' num2str(n) ', k=' num2str(k) ', R=' num2str(k/n,'%.2f') ', ' num2str(ebn0(ii)) ' dB, Dec=' num2str(num_decoded(ii)) ', Err=' num2str(num_errs(ii)) ', Bit_Err=' num2str(num_bit_errs(ii)) ', BLER=' num2str(num_errs(ii)/num_decoded(ii)) ', BER=' num2str(num_bit_errs(ii)/(num_decoded(ii)*n)) ', EG=' num2str(num_queries(ii)/num_decoded(ii),'%.0f')]); 
            % Format and write the text line to the file
            % fprintf(fileID, ...
            % '(end) mlen_bits = 40*k-1 : n=%d, k=%d, R=%.2f, %.1f dB, Dec=%d, Err=%d, Bit_Err=%d, BLER=%.4f, BER=%.6f, EG=%.0f\n', ...
            % n, k, k/n, ebn0(ii), num_decoded(ii), num_errs(ii), num_bit_errs(ii), ...
            % num_errs(ii)/num_decoded(ii), num_bit_errs(ii)/(num_decoded(ii)*n), ...
            % num_queries(ii)/num_decoded(ii));
        end
    end

    % =========================
    % finish statistics for ii
    % =========================

    if isequal(code_class,'ROCB') || isequal(code_class,'ROCB_CRC')
        P_correct(ii)   = num_correct_accepted(ii) / num_messages(ii);
        P_reject(ii)    = num_rejected(ii) / num_messages(ii);
        P_incorrect(ii) = num_incorrect_accepted(ii) / num_messages(ii);
        assert(abs(P_correct(ii) + P_reject(ii) + P_incorrect(ii) - 1) < 1e-12, ...
    'Message-level probabilities do not sum to 1.');
    end

    fprintf(fileID, ...
    '(end) mlen_bits = rand : n=%d, k=%d, R=%.2f, %.1f dB, Dec=%d, Err=%d, Bit_Err=%d, BLER=%.4f, BER=%.6f, EG=%.0f\n', ...
    n, k, k/n, ebn0(ii), num_decoded(ii), num_errs(ii), num_bit_errs(ii), ...
    num_errs(ii)/num_decoded(ii), num_bit_errs(ii)/(num_decoded(ii)*n), ...
    num_queries(ii)/num_decoded(ii));

    % % =========================
    % % SAVE THIS SNR NOW
    % % =========================
    % 
    % if isequal(modulation,'BSC')
    %     filename = ['RESULTS/BSC_' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_' num2str(nmodbits) '.mat'];
    % else
    %     filename = ['RESULTS/BPSK_' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_' num2str(nmodbits) '.mat'];
    % end
    % 
    % % Store the data in a way where results from future simulations with the 
    % % same code can be appended.
    % 
    % if exist(filename, 'file') == 2
    %     S = load(filename,'code');
    %     code = S.code;
    % 
    %     % these = find(code.snr==snr_db(ii));
    %     these = find(abs(code.snr - snr_db(ii)) < 1e-10);
    % 
    %     if isempty(these)
    %         % add new SNR point
    %         code.snr(end+1) = snr_db(ii);
    %         code.num_decoded(end+1)=num_decoded(ii);
    %         code.num_demod_errs(end+1)=num_demod_errs(ii);
    %         code.num_demod_bit_errs(end+1)=num_demod_bit_errs(ii);
    %         code.num_errs(end+1)=num_errs(ii);
    %         code.num_bit_errs(end+1)=num_bit_errs(ii);
    %         code.num_aband(end+1)=num_aband(ii);
    %         code.num_queries(end+1)=num_queries(ii);
    % 
    %         if isequal(code_class,'ROCB') || isequal(code_class,'ROCB_CRC')
    %             code.num_messages(end+1) = num_messages(ii);
    %             code.num_rejected(end+1) = num_rejected(ii);
    %             code.num_correct_accepted(end+1) = num_correct_accepted(ii);
    %             code.num_incorrect_accepted(end+1) = num_incorrect_accepted(ii);
    % 
    %             code.P_correct(end+1) = num_correct_accepted(ii) ./ num_messages(ii);
    %             code.P_reject(end+1) = num_rejected(ii) ./ num_messages(ii);
    %             code.P_incorrect(end+1) = num_incorrect_accepted(ii) ./ num_messages(ii);
    %         end
    %     else
    %         % overwrite existing SNR point
    %         % code.num_decoded(these)=code.num_decoded(these)+num_decoded(ii);
    %         % code.num_demod_errs(these)=code.num_demod_errs(these)+num_demod_errs(ii);
    %         % code.num_demod_bit_errs(these)=code.num_demod_bit_errs(these)+num_demod_bit_errs(ii);
    %         % code.num_errs(these)=code.num_errs(these)+num_errs(ii);
    %         % code.num_bit_errs(these)=code.num_bit_errs(these)+num_bit_errs(ii);
    %         % code.num_aband(these)=code.num_aband(these)+num_aband(ii);
    %         % code.num_queries(these)=code.num_queries(these)+num_queries(ii);
    %         code.num_decoded(these)=num_decoded(ii);
    %         code.num_demod_errs(these)=num_demod_errs(ii);
    %         code.num_demod_bit_errs(these)=num_demod_bit_errs(ii);
    %         code.num_errs(these)=num_errs(ii);
    %         code.num_bit_errs(these)=num_bit_errs(ii);
    %         code.num_aband(these)=num_aband(ii);
    %         code.num_queries(these)=num_queries(ii);
    % 
    %         if isequal(code_class,'ROCB') || isequal(code_class,'ROCB_CRC')
    %             code.num_messages(these) = num_messages(ii);
    %             code.num_rejected(these) = num_rejected(ii);
    %             code.num_correct_accepted(these) = num_correct_accepted(ii);
    %             code.num_incorrect_accepted(these) = num_incorrect_accepted(ii);
    % 
    %             code.P_correct(these) = num_correct_accepted(ii) ./ num_messages(ii);
    %             code.P_reject(these) = num_rejected(ii) ./ num_messages(ii);
    %             code.P_incorrect(these) = num_incorrect_accepted(ii) ./ num_messages(ii);
    %         end
    %     end 
    % else
    %     % First completed SNR point written to this file
    % 
    %     code = struct();
    % 
    %     % =========================
    %     % update code with current SNR
    %     % =========================
    % 
    %     code.class = code_class;
    %     code.n = n;
    %     code.k = k;
    %     code.modulation = modulation;
    %     code.nmodbits = nmodbits;
    % 
    %     code.snr = snr_db(ii);
    % 
    %     code.num_decoded = num_decoded(ii);
    %     code.num_demod_errs = num_demod_errs(ii);
    %     code.num_demod_bit_errs = num_demod_bit_errs(ii);
    %     code.num_errs = num_errs(ii);
    %     code.num_bit_errs = num_bit_errs(ii);
    %     code.num_aband = num_aband(ii);
    %     code.num_queries = num_queries(ii);
    % 
    %     if isequal(code_class,'ROCB') || isequal(code_class,'ROCB_CRC')
    %         code.num_messages = num_messages(ii);
    %         code.num_rejected = num_rejected(ii);
    %         code.num_correct_accepted = num_correct_accepted(ii);
    %         code.num_incorrect_accepted = num_incorrect_accepted(ii);
    % 
    %         code.P_correct = P_correct(ii);
    %         code.P_reject = P_reject(ii);
    %         code.P_incorrect = P_incorrect(ii);
    %     end
    % end
    % 
    % % % =========================
    % % % sort accumulated results
    % % % =========================
    % % 
    % % [~,ord]=sort(code.snr);
    % % 
    % % code.snr = code.snr(ord);
    % % code.num_decoded = code.num_decoded(ord);
    % % code.num_demod_errs = code.num_demod_errs(ord);
    % % code.num_demod_bit_errs = code.num_demod_bit_errs(ord);
    % % code.num_errs = code.num_errs(ord);
    % % code.num_bit_errs = code.num_bit_errs(ord);
    % % code.num_aband = code.num_aband(ord);
    % % code.num_queries = code.num_queries(ord);
    % % if isequal(code_class,'ROCB') || isequal(code_class,'ROCB_CRC')
    % %     code.num_messages = code.num_messages(ord);
    % %     code.num_rejected = code.num_rejected(ord);
    % %     code.num_correct_accepted = code.num_correct_accepted(ord);
    % %     code.num_incorrect_accepted = code.num_incorrect_accepted(ord);
    % % 
    % %     code.P_correct = code.P_correct(ord);
    % %     code.P_reject = code.P_reject(ord);
    % %     code.P_incorrect = code.P_incorrect(ord);
    % % end
    % 
    % 
    % 
    % % =========================
    % % recompute summary values
    % % =========================
    % 
    % code.ebn0 = code.snr-10*log10(k/n)-10*log10(nmodbits); %Eb/N0
    % code.R = code.k/code.n; % Code rate
    % code.BLERdemod = code.num_demod_errs./code.num_decoded; % Demod BLER
    % code.BERdemod = code.num_demod_bit_errs./(code.num_decoded*code.n); % Demod BER
    % code.BLER = code.num_errs./code.num_decoded; % Decoded BLER
    % code.BER = code.num_bit_errs./(code.num_decoded*code.n); % Decoded BER
    % code.EG = code.num_queries./code.num_decoded; % Average number of code-book queries
    % code.max_query = max_query;
    % 
    % 
    % % =========================
    % % write immediately
    % % =========================
    % 
    % save(filename,'code')
    % 
    % if isequal(modulation,'BSC')
    %     disp(['BSC; Decoder=' DECODER ', code=' code.class ', ' num2str(sum(num_decoded)) ' packets decoded, ebn0(ii)= ' mat2str(ebn0(ii)) '.'])
    % 
    %     fprintf(fileID, ...
    %         'BSC Decoder=%s, code=%s, %d packets decoded.', ...
    %         DECODER, code.class, sum(num_decoded));
    % else
    %     disp(['Decoder=' DECODER ', code=' code.class ', ' num2str(sum(num_decoded)) ' packets decoded, ebn0(ii)= ' mat2str(ebn0(ii)) '.'])
    % 
    %     fprintf(fileID, ...
    %         'Decoder=%s, code=%s, %d packets decoded.', ...
    %         DECODER, code.class, sum(num_decoded));
    % end

end


% =========================
% SAVE THE RESULTS
% =========================

addpath('RESULTS/');

if isequal(modulation,'BSC')
    filename = ['RESULTS/BSC_' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_' num2str(nmodbits) '.mat'];
else
    filename = ['RESULTS/BPSK_' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_' num2str(nmodbits) '.mat'];
end

% Store the data in a way where results from future simulations with the 
% same code can be appended.

if exist(filename, 'file') == 2
    load(filename,'code')

    % Backward compatibility with older result files - Older result files may not contain message-level statistics
    if isequal(code_class,'ROCB') || isequal(code_class,'ROCB_CRC')
        if ~isfield(code,'num_messages')
            code.num_messages = nan(size(code.snr));
        end
        if ~isfield(code,'num_rejected')
            code.num_rejected = nan(size(code.snr));
        end
        if ~isfield(code,'num_correct_accepted')
            code.num_correct_accepted = nan(size(code.snr));
        end
        if ~isfield(code,'num_incorrect_accepted')
            code.num_incorrect_accepted = nan(size(code.snr));
        end
    end

    for ii=1:length(snr_db)
        these = find(code.snr==snr_db(ii));
        if isempty(these)
            code.snr(end+1) = snr_db(ii);
            code.num_decoded(end+1)=num_decoded(ii);
            code.num_demod_errs(end+1)=num_demod_errs(ii);
            code.num_demod_bit_errs(end+1)=num_demod_bit_errs(ii);
            code.num_errs(end+1)=num_errs(ii);
            code.num_bit_errs(end+1)=num_bit_errs(ii);
            code.num_aband(end+1)=num_aband(ii);
            code.num_queries(end+1)=num_queries(ii);
            if isequal(code_class,'ROCB') || isequal(code_class,'ROCB_CRC')
                code.num_messages(end+1)           = num_messages(ii);
                code.num_rejected(end+1)           = num_rejected(ii);
                code.num_correct_accepted(end+1)   = num_correct_accepted(ii);
                code.num_incorrect_accepted(end+1) = num_incorrect_accepted(ii);
            end
        else
            % code.num_decoded(these)=code.num_decoded(these)+num_decoded(ii);
            % code.num_demod_errs(these)=code.num_demod_errs(these)+num_demod_errs(ii);
            % code.num_demod_bit_errs(these)=code.num_demod_bit_errs(these)+num_demod_bit_errs(ii);
            % code.num_errs(these)=code.num_errs(these)+num_errs(ii);
            % code.num_bit_errs(these)=code.num_bit_errs(these)+num_bit_errs(ii);
            % code.num_aband(these)=code.num_aband(these)+num_aband(ii);
            % code.num_queries(these)=code.num_queries(these)+num_queries(ii);
            code.num_decoded(these)=num_decoded(ii);
            code.num_demod_errs(these)=num_demod_errs(ii);
            code.num_demod_bit_errs(these)=num_demod_bit_errs(ii);
            code.num_errs(these)=num_errs(ii);
            code.num_bit_errs(these)=num_bit_errs(ii);
            code.num_aband(these)=num_aband(ii);
            code.num_queries(these)=num_queries(ii);
            if isequal(code_class,'ROCB') || isequal(code_class,'ROCB_CRC')
                code.num_messages(these)           = num_messages(ii);
                code.num_rejected(these)           = num_rejected(ii);
                code.num_correct_accepted(these)   = num_correct_accepted(ii);
                code.num_incorrect_accepted(these) = num_incorrect_accepted(ii);
            end
        end
    end     
else
    % =========================
    % update code with current SNR
    % =========================

    code.n = n;
    code.k = k;
    code.modulation = modulation;
    code.nmodbits = nmodbits;
    code.snr = snr_db;
    code.num_decoded=num_decoded;
    code.num_demod_errs=num_demod_errs;
    code.num_demod_bit_errs=num_demod_bit_errs;
    code.num_errs=num_errs;
    code.num_bit_errs=num_bit_errs;
    code.num_aband=num_aband;
    code.num_queries=num_queries;
    if isequal(code_class,'ROCB') || isequal(code_class,'ROCB_CRC')
        code.num_messages           = num_messages;
        code.num_rejected           = num_rejected;
        code.num_correct_accepted   = num_correct_accepted;
        code.num_incorrect_accepted = num_incorrect_accepted;
    end
end

% =========================
% sort accumulated results
% =========================

% Sort according to increasing snr;
[~,ord]=sort(code.snr);

code.snr = code.snr(ord);
code.num_decoded = code.num_decoded(ord);
code.num_demod_errs = code.num_demod_errs(ord);
code.num_demod_bit_errs = code.num_demod_bit_errs(ord);
code.num_errs = code.num_errs(ord);
code.num_bit_errs = code.num_bit_errs(ord);
code.num_aband = code.num_aband(ord);
code.num_queries = code.num_queries(ord);
if isequal(code_class,'ROCB') || isequal(code_class,'ROCB_CRC')
    if ~isfield(code,'num_messages')
        code.num_messages = nan(size(code.snr));
    end
    if ~isfield(code,'num_rejected')
        code.num_rejected = nan(size(code.snr));
    end
    if ~isfield(code,'num_correct_accepted')
        code.num_correct_accepted = nan(size(code.snr));
    end
    if ~isfield(code,'num_incorrect_accepted')
        code.num_incorrect_accepted = nan(size(code.snr));
    end

    code.num_messages           = code.num_messages(ord);
    code.num_rejected           = code.num_rejected(ord);
    code.num_correct_accepted   = code.num_correct_accepted(ord);
    code.num_incorrect_accepted = code.num_incorrect_accepted(ord);
end

% =========================
% recompute summary values
% =========================
% Summary statistics   
code.ebn0 = code.snr-10*log10(k/n)-10*log10(nmodbits); %Eb/N0
code.R = code.k/code.n; % Code rate
code.BLERdemod = code.num_demod_errs./code.num_decoded; % Demod BLER
code.BERdemod = code.num_demod_bit_errs./(code.num_decoded*code.n); % Demod BER
code.BLER = code.num_errs./code.num_decoded; % Decoded BLER
code.BER = code.num_bit_errs./(code.num_decoded*code.n); % Decoded BER
code.EG = code.num_queries./code.num_decoded; % Average number of code-book queries
code.max_query = max_query;

if isequal(code_class,'ROCB') || isequal(code_class,'ROCB_CRC')
    code.P_correct   = code.num_correct_accepted   ./ code.num_messages;
    code.P_reject    = code.num_rejected           ./ code.num_messages;
    code.P_incorrect = code.num_incorrect_accepted ./ code.num_messages;

    valid_msg_stats = ...
        ~isnan(code.num_messages) & ...
        ~isnan(code.num_rejected) & ...
        ~isnan(code.num_correct_accepted) & ...
        ~isnan(code.num_incorrect_accepted);

    assert(all(abs( ...
        code.P_correct(valid_msg_stats) + ...
        code.P_reject(valid_msg_stats) + ...
        code.P_incorrect(valid_msg_stats) - 1) < 1e-12), ...
        'Saved message-level probabilities do not sum to 1.');
end

% =========================
% write immediately
% =========================
save(filename,'code')

if isequal(modulation,'BSC')
    disp(['BSC; Decoder=' DECODER ', code=' code.class ', ' num2str(sum(num_decoded)) ' packets decoded.'])

    fprintf(fileID, ...
        'BSC Decoder=%s, code=%s, %d packets decoded.', ...
        DECODER, code.class, sum(num_decoded));
else
    disp(['Decoder=' DECODER ', code=' code.class ', ' num2str(sum(num_decoded)) ' packets decoded.'])

    fprintf(fileID, ...
        'Decoder=%s, code=%s, %d packets decoded.', ...
        DECODER, code.class, sum(num_decoded));
end

% Only close Report.txt after ALL SNR points are complete
% Close the file
fclose(fileID);



