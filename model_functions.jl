function get_matches(Data, f23_df)
##################################
    
    # Create a dictionary to store indices of hex strings in Data
    index_dict = Dict{String, Vector{Int}}()
    
    # Populate the dictionary
    for (i, hex_str) in enumerate(Data)
        if haskey(index_dict, hex_str)
            push!(index_dict[hex_str], i)
        else
            index_dict[hex_str] = [i]
        end
    end
    
    # Initialize a vector to store indices
    matching_indices = []
    
    # Iterate through each hex string in f23_df and lookup in the dictionary
    for hex_str in f23_df.Match_vector
        if haskey(index_dict, hex_str)
            push!(matching_indices, index_dict[hex_str][1])
        else
            push!(matching_indices, nothing)  # If no match found, store an empty vector
        end
    end

    f23_df[!,"Data_vector"] = matching_indices

    return(f23_df)

end    # get_matches()


# function to calculate selected parameters from Spectrum synchronisation message (0xF23)
function process_f23(f23_vals)
##############################
    
    # refer to DWTP (Ver. 16 January2019) Section 4.3 pp.43-44

    # get Timestamp in UTC
    timestamp = unix2datetime(parse(Int, bitstring(f23_vals[3]) * bitstring(f23_vals[4]) * bitstring(f23_vals[5]) * bitstring(f23_vals[6]); base=2))
    
    # convert time to Australian Eastern Standard Time
    timestamp = timestamp + Hour(0)  # Adjust this for the correct time zone

    # get Data Stamp
    data_stamp = parse(Int, bitstring(f23_vals[7]) * bitstring(f23_vals[8]); base=2)

    # get Segments Used
    segments_used = parse(Int, bitstring(f23_vals[9]) * bitstring(f23_vals[10]) * bitstring(f23_vals[11]); base=2)

    # get Sample Number
    sample_number = parse(Int, bitstring(f23_vals[12]) * bitstring(f23_vals[13]); base=2)

    # Create Match Vector
    match_vector = lpad(string(f23_vals[14], base=16), 2, "0")
    for i in 15:22
        match_vector = match_vector * lpad(string(f23_vals[i], base=16), 2, "0")
    end
    
    return(timestamp, segments_used, match_vector, sample_number)
    
end    #  process_f23()


# convert binary data into F23_df and Hex array
function get_hex_array(infil)
#############################
    
    # Read binary data from the input file
    println("Reading BINARY data from ", infil)
    flush(stdout)
    data = reinterpret(UInt8, read(infil))
    
    # Turn the data vector into a matrix of 12 values matching hexadecimal bytes
    cols = 12
    rows = Int(length(data) / cols)
    mat = reshape(view(data, :), cols, :)
    
    # Calculate the Heave, North, and West displacements
    hex_matrix = string.(mat'[:,1:9], base=16, pad=2)
    Data = [join(row) for row in eachrow(hex_matrix)]
    
    println("All file data read!")
    
    # Interleave the last 3 matrix columns (10, 11, 12) to form the packet vector
    packet = collect(Iterators.flatten(zip(mat[10,:], mat[11,:], mat[12,:])))
    
    # Find all occurrences of 0x7e in the packet vector
    aa = findall(x -> x == 0x7e, vec(packet))
    
    # Create DataFrame to hold the processed data
    f23_df = DataFrame(Date = DateTime[], Segments = Int[], Match_vector = String[], Sample_number = Int[])
    
    # Decode the packet data into messages
    max_val = length(aa) - 1
    
    for i in 1:max_val
        first = aa[i] + 1
        last = aa[i + 1]
        
        if (last - first > 1)
            decoded = packet[first:last-1]
            
            # Handle the 0x7d escape sequences (XOR with 0x20)
            bb = findall(x -> x == 0x7d, decoded)
            for ii in bb
                decoded[ii + 1] = decoded[ii + 1] ⊻ 0x20
            end
            deleteat!(decoded, bb)
            
            # If the message is F23 (0x23)
            if decoded[2] == 0x23
                timestamp, segments_used, match_vector, sample_number = process_f23(decoded)
                push!(f23_df, [timestamp, segments_used, match_vector, sample_number])
            end
        end
    end
    
    # Remove duplicates from f23_df
    f23_df = unique(f23_df);

    return(f23_df, Data)
    
end    # get_hex_array()


function get_start_end_dates(f23_df,found_list)
###############################################
    
    start_date = f23_df[found_list[1],:].Date - Minute(30) # <------- NOTE subtracted 30min from start_date to match Waves4 results
    segments = f23_df[found_list[1],:].Segments
#   match_vector = f23_df[found_list[1],:].Match_vector
    sample_nos = f23_df[found_list[1],:].Sample_number
    data_vector = f23_df[found_list[1],:].Data_vector
    start_val = data_vector - Int(sample_nos/2) + 1
    end_val = data_vector
    
    return(start_date,start_val, end_val)
    
end    #(get_start_end_dates)
    
  
function get_displacement(Data, start_val, end_val)
###################################################
# Decode the real time data to displacements - See DWTP (16 Jan 2019) 2.1.1 p. 19    
    
    arry = collect(Iterators.flatten(zip(SubString.(Data, start_val, end_val),SubString.(Data, start_val+9, end_val+9))));
    displacements = [parse(Int, SubString.(i, 1, 1), base=16)*16^2 + parse(Int, SubString.(i, 2, 2), base=16)*16^1 + parse(Int, SubString.(i, 3, 3), base=16)*16^0 for i in arry]    
    
    displacements[findall(>=(2048), displacements)] = displacements[findall(>=(2048), displacements)] .- 4096;
    displacements = 0.457*sinh.(displacements/457)    # see DWTP p.19 (16)
   
    return(displacements)
    
end    # get_displacement()


function get_hnw(Data,start_val,end_val)
######################################## 
    
    # get WSEs for desired 30-minute record
    heave = get_displacement(Data[start_val:end_val,:], 1, 3);              
    north = get_displacement(Data[start_val:end_val,:], 4, 6);
    west = get_displacement(Data[start_val:end_val,:], 7, 9);
    
    # Check for missing or extra points in data
    for wse in [heave, north, west]
        
        wse_length = length(wse)
        
        if wse_length > REC_LENGTH

            # truncate if too long
            wse = wse[1:REC_LENGTH]
            
        else

            # zero pad if too short (leave it unchanged if right length)
            append!(wse,zeros(REC_LENGTH-wse_length))
            
        end      

    end
    
    return (heave, north, west)
    
end    # get_hnw()


# Function to calculate confidence limits
function calc_confidence_limits(data, confidence_interval)
##########################################################
    
    mean_val = mean(data)
    std_dev = std(data)
    upper_limit = mean_val + confidence_interval * std_dev
    lower_limit = mean_val - confidence_interval * std_dev
    
    return (lower_limit, upper_limit)
    
end    # calc_confidence_limits()


# Function to compute modified z-scores and find outliers
function modified_z_score(data, threshold)
##########################################
    
    med = median(data)
    mad = median(abs.(data .- med))
    mod_z_scores = 0.6745 * (data .- med) ./ mad
    outlier_indices = findall(x -> abs(x) > threshold, mod_z_scores)

    return(outlier_indices, mod_z_scores)
    
end    # modified_z_score()


# Function for dynamic threshold based on mean wave height
function RETIRED_dynamic_z_score_threshold(heave, base_threshold=3.0, k=0.5)
####################################################################
#==
    mean_wave_height = mean(heave)
    std_wave_height = std(heave)
    dynamic_threshold = base_threshold * (1 + k * (mean_wave_height / std_wave_height))
==#
    max_threshold = 3.29
    
    heave_range = maximum(heave) - minimum(heave)
    std_wave_height = std(heave)

    # Scale the threshold based on the range of heave values
    # Adjust the scaling factor (e.g., 0.1) as needed to fit your data
    scaling_factor = 0.1 * (heave_range / std_wave_height)

    # Calculate the dynamic threshold
    dynamic_threshold = base_threshold + scaling_factor

    # Clamp the threshold between the defined limits
    dynamic_threshold = clamp(dynamic_threshold, base_threshold, max_threshold)
    
    return(dynamic_threshold)
    
end    # dynamic_z_score_threshold()


# Function for dynamic threshold based on wave amplitudes
function dynamic_z_score_threshold(heave, base_threshold=3.0, k=0.5)
####################################################################
    
    # Calculate the amplitude of heave
    amplitude = abs.(heave)
    
    # Calculate mean and standard deviation of the amplitudes
    mean_amplitude = mean(amplitude)
    std_amplitude = std(amplitude)
    
    # Calculate the range of amplitudes
    range_amplitude = maximum(amplitude) - minimum(amplitude)
    
    # Calculate dynamic threshold
    dynamic_threshold = base_threshold * (1 + k * (range_amplitude / std_amplitude))
    
    # Cap the dynamic threshold
    ############################################################################################################
    ##    confidence_interval = 2.576  # corresponds to a 99% confidence interval (for a normal distribution)
    ##    confidence_interval = 3.0    # corresponds to a 99.73% confidence interval (for a normal distribution)    
    ##    confidence_interval = 3.29   # corresponds to a 99.9% confidence interval (for a normal distribution)
    ############################################################################################################
    max_threshold = 3.29  # see above
    dynamic_threshold = min(dynamic_threshold, max_threshold)
    
    return dynamic_threshold
end  # dynamic_z_score_threshold()


function pad_or_truncate(record, target_length=REC_LENGTH)
##########################################################

    length(record) < target_length ? vcat(record, zeros(Float32, target_length - length(record))) :
                                     record[1:target_length]

end    # pad_or_truncate()


function get_heave_north_west(Data, f23_df)
###########################################
    
    REC_LENGTH = 4608
    num_records = nrow(f23_df)

    # Initialize a 3D matrix for heave, north, and west data, and a vector for datetime values
    hnw_array = Array{Float32}(undef, REC_LENGTH, num_records, 3)
    X_date = Vector{DateTime}(undef, num_records)

    println("Calculating Heave, North, and West values now!")
    
    record_count = 0  # Counter for valid records

    for idx in 1:num_records
        if !isnothing(f23_df.Data_vector[idx])
            start_date, start_val, end_val = get_start_end_dates(f23_df, idx)

            if start_val > 0
                print(".")
                
                # Get heave, north, and west values for this record
                heave, north, west = get_hnw(Data, start_val, end_val)
                
                # Pad or truncate to ensure we have REC_LENGTH points
                hnw_array[:, record_count + 1, 1] = pad_or_truncate(heave, REC_LENGTH)
                hnw_array[:, record_count + 1, 2] = pad_or_truncate(north, REC_LENGTH)
                hnw_array[:, record_count + 1, 3] = pad_or_truncate(west, REC_LENGTH)
                
                # Store the start date
                X_date[record_count + 1] = start_date
                record_count += 1
            end

        end
    end

    # Resize X_date to only contain valid entries
    X_date = X_date[1:record_count]
    hnw_array = hnw_array[:, 1:record_count, :]  # Resize to remove unused columns

    return(hnw_array, X_date)

end    #  get_heave()


function get_heave(Data, f23_df)
################################
    
    heave_array = []
    X_date = []

    println("Calculating Heave values now!")
    
    for idx in 1:nrow(f23_df)

        if !isnothing(f23_df.Data_vector[idx])
    
            start_date, start_val, end_val = get_start_end_dates(f23_df,idx)
            if start_val > 0
                print(".")
                heave, north, west = get_hnw(Data,start_val,end_val)

                # ensure we have REC_LENGTH data points
                push!(heave_array,pad_or_truncate(heave, REC_LENGTH))
                push!(X_date,start_date)
            end

        end
    
    end

    return(hcat(heave_array...), X_date)

end    # get_heave()


# Need to check first row of the f23_df in case 23:00 is stored there
function f23_first_row_check(f23_df)
####################################
    
    # Get the first row of the DataFrame
    first_row = first(f23_df)
    
    # Check if the time of the first row's Date column is 23:00:00
    time_of_first_row = Time(first_row.Date)

    if time_of_first_row == Time(23, 0, 0)

        if ismissing(first_row.Data_vector) || isnothing(first_row.Data_vector) || isnan(first_row.Data_vector)
            f23_df = f23_df[2:end, :]  # Drop the first row
        end

    end
    
    return(f23_df)
    
end    # f23_first_row_check()


function min_max_normalize_matrix(X)
####################################
    
    min_vals = minimum(X, dims=1)  # Compute min for each column
    max_vals = maximum(X, dims=1)  # Compute max for each column
    
    return((X .- min_vals) ./ (max_vals .- min_vals))
    
end    # min_max_normalize_matrix()


function calc_reconstruction_errors(data_matrix, model)
#######################################################
    
    reconstruction_errors = Float32[]
    
    for record in eachcol(data_matrix)
        reconstructed_record = model(record)
        error = mean((reconstructed_record .- record).^2)
        push!(reconstruction_errors, error)
    end
    
    return(reconstruction_errors)
    
end    # calc_reconstruction_errors()


function detect_outliers(X_data, X_date, training_data_good, training_data_bad, model, weight_factor=2.0)
#########################################################################################################
    
    # Normalize training data
    errors_good = calc_reconstruction_errors(Float32.(min_max_normalize_matrix(training_data_good)), model)
    errors_bad = calc_reconstruction_errors(Float32.(min_max_normalize_matrix(training_data_bad)), model)

    # Calculate normalization parameters
    mean_train_selected = mean(training_data_good, dims=1)[:, 1:size(X_data)[2]]
    std_train_selected = std(training_data_good, dims=1)[:, 1:size(X_data)[2]]
    
    # Normalize new data
    X_data_normalized = min_max_normalize_matrix(X_data)

    # Predict and calculate reconstruction error
    predicted_X_data = model(X_data_normalized)
    reconstruction_error = sum((X_data_normalized .- predicted_X_data) .^ 2, dims=1)
    reconstruction_error_vector = vec(reconstruction_error)
    normalized_reconstruction_error = min_max_normalize_matrix(reconstruction_error_vector)
    inverted_reconstruction_error = 1.0 .- normalized_reconstruction_error

    # Adaptive thresholding based on the good and bad data thresholds
    weighted_errors_bad = errors_bad .* weight_factor
    good_threshold = quantile(errors_good, 0.95)
    bad_threshold = quantile(weighted_errors_bad, 0.995)
    
    # Adaptive uncertain and bad thresholds using median and std
    inverted_median = median(inverted_reconstruction_error)
    inverted_std = std(inverted_reconstruction_error)
    uncertain_threshold = inverted_median + inverted_std
    bad_threshold = inverted_median + 1.5 * inverted_std

    # Identify outliers and uncertain points based on thresholds
    outliers = findall(inverted_reconstruction_error .> bad_threshold)
    uncertain_indices = findall(x -> uncertain_threshold < x <= bad_threshold, inverted_reconstruction_error)

    # Map indices to dates
    outlier_dates = X_date[outliers]
    uncertain_dates = X_date[uncertain_indices]

    return(outliers, uncertain_indices, outlier_dates, uncertain_dates, good_threshold, bad_threshold)

end    # detect_outliers()


function do_plots(ii, start_time, heave)
########################################
    
    end_time = start_time + Minute(30)
    xvals = start_time + Microsecond.((0:REC_LENGTH-1) / SAMPLE_FREQUENCY * 1000000)

    Q1 = quantile(heave, 0.25)
    Q3 = quantile(heave, 0.75)

    multiplier = 1.5

    IQR = Q3 - Q1
    lower_bound = Q1 - multiplier * IQR
    upper_bound = Q3 + multiplier * IQR

    # Plot initialization
    p1 = plot(size=(2000, 400), dpi=100, framestyle=:box, fg_legend=:transparent, bg_legend=:transparent, 
        legend=:topright, xtickfont=font(8), ytickfont=font(8), bottommargin=5Plots.mm, 
        grid=true, gridlinewidth=0.125, gridstyle=:dot, gridcolor=:grey, gridalpha=0.5)
    
    p1 = hspan!([lower_bound, upper_bound], fillcolor=:lightblue, fillalpha=:0.125, label="IQR limits")
    
    tm_tick = range(start_time, end_time, step=Minute(1))
    ticks = Dates.format.(tm_tick, "MM")
    
    # Calculate dynamic confidence interval
    confidence_interval = 3.29 # threshold at the 99.9th percentile level

    # Identify z_scores using modified z-score
    z_score_indices, mod_z_scores = modified_z_score(heave, confidence_interval)
    if !isempty(z_score_indices)
        scatter!(p1, xvals[z_score_indices], heave[z_score_indices], 
            markersize=4, markerstrokecolor=:red, markerstrokewidth=1, 
            markercolor=:white, markershape=:circle, label="Modified Z-score beyond 99.9% confidence limits")
    end

    # Plot confidence limits
    confidence_limits = calc_confidence_limits(heave, confidence_interval)
    hline!(p1, [confidence_limits[1], confidence_limits[2]], color=:red, lw=0.5, linestyle=:dash, label="99.9% confidence limits")

    # Plot heave data
    plot!(p1, xvals, heave, xlims=(xvals[1], xvals[end]), lw=0.5, lc=:blue, alpha=0.5, 
        xticks=(tm_tick, ticks), label="")

    # Annotate plot with the number of outliers and confidence interval
    num_outliers = length(z_score_indices)
    suspect_string = string("  ", string(ii)," ",Dates.format(start_time, "yyyy-mm-dd HH:MM"), " - ", num_outliers, " Possible outliers") # using Confidence Interval of ", 
##        @sprintf("%.2f", confidence_interval))
    annotate!(p1, xvals[1], maximum(heave) * 0.9, text(suspect_string, :left, 10, :blue))

    display(p1)
    
end    # do_plots()  


function calc_and_plot_bounds(xvals, heave, lc, hspan_color, label)
###################################################################

    Q1 = quantile(heave, 0.25)
    Q3 = quantile(heave, 0.75)

    multiplier = 1.5

    IQR = Q3 - Q1
    lower_bound = Q1 - multiplier * IQR
    upper_bound = Q3 + multiplier * IQR

    # Calculate dynamic confidence interval
    confidence_interval = 3.29 # threshold at the 99.9th percentile level

    # Identify z_scores using modified z-score
    z_score_indices, mod_z_scores = modified_z_score(heave, confidence_interval)

    # Plot confidence limits
    confidence_limits = calc_confidence_limits(heave, confidence_interval)

    tm_tick = range(xvals[1], xvals[end], step=Minute(1))
    ticks = Dates.format.(tm_tick, "MM")


    px = plot(xvals, heave, xlims=(xvals[1], xvals[end]), lw=0.5, lc=lc, alpha=0.5, 
        xticks=(tm_tick, ticks), label=label, legendfontsize=12)

    if !isempty(z_score_indices)
        scatter!(px, xvals[z_score_indices], heave[z_score_indices], 
            markersize=4, markerstrokecolor=:red, markerstrokewidth=1, 
            markercolor=:white, markershape=:circle, label="Modified Z-score beyond 99.9% confidence limits")
    end

    px = hspan!([lower_bound, upper_bound], fillcolor=hspan_color, fillalpha=:0.125, label="IQR limits")
    px = hline!([confidence_limits[1], confidence_limits[2]], color=:red, lw=0.5, linestyle=:dash, label="99.9% confidence limits")            
    
    return(px)

end    # calc_and_plot_bounds()


function do_heave_north_west_plots(ii, start_time, X_data)
##########################################################

    heave = X_data[:, ii, 1]
    north = X_data[:, ii, 2]
    west = X_data[:, ii, 3]
    
    end_time = start_time + Minute(30)
    xvals = start_time + Microsecond.((0:REC_LENGTH-1) / SAMPLE_FREQUENCY * 1000000)
   
    tm_tick = range(start_time, end_time, step=Minute(1))
    ticks = Dates.format.(tm_tick, "MM")
            
    p1 = calc_and_plot_bounds(xvals, heave, :blue, :lightblue, "Heave")

    p2 = calc_and_plot_bounds(xvals, north, :red, :pink, "North")

    p3 = calc_and_plot_bounds(xvals, west, :green, :lightgreen, "West")

    # Annotate plot with the number of outliers and confidence interval
##    num_outliers = length(z_score_indices)
##    suspect_string = string("  ", string(ii)," ",Dates.format(start_time, "yyyy-mm-dd HH:MM"), " - ", num_outliers, " Possible outliers") # using Confidence Interval of ", 
##        @sprintf("%.2f", confidence_interval))
##    annotate!(p1, xvals[1], maximum(heave) * 0.9, text(suspect_string, :left, 10, :blue))

    date_string = Dates.format(start_time, "yyyy-mm-dd HH:MM") * " UTC"    # <- NOTE: .bva times are UTC!

    plot_displacements = plot(p1, p2, p3, size=(2000, 1000), layout=(3,1), dpi=100, framestyle=:box, fg_legend=:transparent, bg_legend=:transparent, 
    legend=:topright, xtickfont=font(8), ytickfont=font(8), bottommargin=5Plots.mm, suptitle=date_string,
    grid=true, gridlinewidth=0.125, gridstyle=:dot, gridcolor=:grey, gridalpha=0.5)
 
    display(plot_displacements)
    
end    # do_heave_north_west_plots()  


