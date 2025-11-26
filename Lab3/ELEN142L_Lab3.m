%% Part 1
clear
clc
snr_db = zeros(1,10);
for i=1:1:10
    n_loops = 0;
    n_symbol_errors = 0;
    n_bit_errors = 0;
    n_block_errors = 0;
    snr_db(i) = 4 + i; % 5 6 7 8 ...
    snr(i) = 10.^(snr_db(i)/10);
    
    while n_symbol_errors < 100
        tx_signal(1:100,1) = sign(randn([100,1]));
        noise = 1/sqrt(snr(i))*randn(100,1);
        
        rx_signal = tx_signal + noise;
        rx_decoded = sign(rx_signal);
        
        cur_n_symbol_errors = sum(tx_signal ~= rx_decoded);
        n_symbol_errors = n_symbol_errors + cur_n_symbol_errors;
        n_bit_errors = n_bit_errors + cur_n_symbol_errors;
        if cur_n_symbol_errors > 0
            n_block_errors = n_block_errors + 1;
        end
        n_loops = n_loops + 1;
    end
    symbol_error_rate(i) = n_symbol_errors / (n_loops * 100);
    bit_error_rate(i) = n_bit_errors / (n_loops * 100);
    block_error_rate(i) = n_block_errors / n_loops;
    disp(i);
end


semilogy(snr_db,bit_error_rate,snr_db,block_error_rate,'--');
legend('Bit Error Rate','Block Error Rate (100 bits)')
xlabel(' SNR (dB)')
ylabel('Error Rate')
grid;

%% Part 2
clear
clc
snr_db = zeros(1,10);
QPSK = [ 1+1j, 1-1j, -1+1j, -1-1j ];
for i=1:1:10
    n_loops = 0;
    n_symbol_errors = 0;
    n_bit_errors = 0;
    n_block_errors = 0;
    snr_db(i) = 4 + i; % 5 6 7 8 ...
    snr(i) = 10.^(snr_db(i)/10);
    
    while n_symbol_errors < 100
        QPSK_indices = randi([1,4], 100, 1);
        tx_signal = QPSK(QPSK_indices);
        noise = 1/sqrt(snr(i))*(randn(1,100) + 1j*randn(1,100));
        
        rx_signal = tx_signal + noise;
        decode_real = sign(real(rx_signal));
        decode_imag = sign(imag(rx_signal));
        
        real_errors = (decode_real ~= real(tx_signal));
        imag_errors = (decode_imag ~= imag(tx_signal));

        n_bit_errors = n_bit_errors + sum(real_errors) + sum(imag_errors);
        cur_n_symbol_errors = sum(real_errors | imag_errors);
        n_symbol_errors = n_symbol_errors + cur_n_symbol_errors;
        n_block_errors = n_block_errors + sign(cur_n_symbol_errors);

        n_loops = n_loops + 1;
    end
    symbol_error_rate(i) = n_symbol_errors / (n_loops * 100);
    bit_error_rate(i) = n_bit_errors / (n_loops * 100 * 2);
    block_error_rate(i) = n_block_errors / n_loops;
    disp(i);
end

semilogy(snr_db,symbol_error_rate,snr_db,bit_error_rate,'--',snr_db,block_error_rate,'-o')
legend('Symbol Error Rate','Bit Error Rate','Block Error Rate (100 bits)')
xlabel(' SNR (dB)')
ylabel('Error Rate')
grid;

%% Part 3
clear
clc
time_vec = linspace(0, 0.005, 100);
doppler = (0.3 * 2 * 10e9 ) / (3 * 10e8);
snr_a_db = zeros(1,10);
snr_a = 0;
QPSK = [ 1+1j, 1-1j, -1+1j, -1-1j ];
for i=1:1:10
    n_loops = 0;
    n_symbol_errors = 0;
    n_bit_errors = 0;
    n_block_errors = 0;
    snr_db(i) = 4 + i; % 5 6 7 8 ...
    snr(i) = 10.^(snr_db(i)/10);
    
    while n_symbol_errors < 10000
        alpha = Jakes(time_vec, 1, doppler);
        snr_a = snr_a + mean(abs(alpha).^2)*snr(i);
        QPSK_indices = randi([1,4], 100, 1);
        tx_signal = QPSK(QPSK_indices);
        noise = 1/sqrt(snr(i))*(randn(1,100) + 1j*randn(1,100));
        
        rx_signal = alpha.'.*tx_signal;
        rx_signal = rx_signal + noise;
        
        decode_real = sign(real(rx_signal./alpha.'));
        decode_imag = sign(imag(rx_signal./alpha.'));
        
        real_errors = (decode_real ~= real(tx_signal));
        imag_errors = (decode_imag ~= imag(tx_signal));

        n_bit_errors = n_bit_errors + sum(real_errors) + sum(imag_errors);
        cur_n_symbol_errors = sum(real_errors | imag_errors);
        n_symbol_errors = n_symbol_errors + cur_n_symbol_errors;
        n_block_errors = n_block_errors + sign(cur_n_symbol_errors);

        n_loops = n_loops + 1;
    end
    symbol_error_rate(i) = n_symbol_errors / (n_loops * 100);
    bit_error_rate(i) = n_bit_errors / (n_loops * 100 * 2);
    block_error_rate(i) = n_block_errors / n_loops;
    snr_a_db(i) = 10*log10(snr_a / n_loops);
    disp(i);
end

semilogy(snr_a_db,symbol_error_rate,snr_a_db,bit_error_rate,'--',snr_a_db,block_error_rate,'-o')
legend('Symbol Error Rate','Bit Error Rate','Block Error Rate (100 bits)')
xlabel(' SNR (dB)')
ylabel('Error Rate')
grid;
