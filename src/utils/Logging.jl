module Logging

export Logger, log_info, log_warn, log_error, log_event, log_states, init_logger, close_logger

using Dates
using ..ComponentIndex

"""
    Logger

Structure to manage logging functionality.
"""
mutable struct Logger
    file_handle::Union{IOStream, Nothing}
    state_file_handle::Union{IOStream, Nothing}
    console_output::Bool
    log_level::Symbol
    
    # Constructor
    function Logger(
        log_file::Union{String, Nothing}=nothing,
        state_log_file::Union{String, Nothing}=nothing;
        console_output::Bool=true,
        log_level::Symbol=:info
    )
        # Create file handles if paths are provided
        file_handle = nothing
        state_file_handle = nothing
        
        if log_file !== nothing
            # Create directory if it doesn't exist
            dir = dirname(log_file)
            if !isdir(dir) && !isempty(dir)
                mkpath(dir)
            end
            file_handle = open(log_file, "w")
        end
        
        if state_log_file !== nothing
            # Create directory if it doesn't exist
            dir = dirname(state_log_file)
            if !isdir(dir) && !isempty(dir)
                mkpath(dir)
            end
            state_file_handle = open(state_log_file, "w")
            
            # Write header for state log file
            write(state_file_handle, "time,variable,value\n")
        end
        
        return new(file_handle, state_file_handle, console_output, log_level)
    end
end

"""
    init_logger(log_dir::String)

Initialize a logger with standard log files in the specified directory.
"""
function init_logger(log_dir::String="logs")
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    
    # Create directory if it doesn't exist
    if !isdir(log_dir)
        mkpath(log_dir)
    end
    
    # Create log files with timestamp
    log_file = joinpath(log_dir, "simulation_$(timestamp).log")
    state_log_file = joinpath(log_dir, "states_$(timestamp).csv")
    
    return Logger(log_file, state_log_file)
end

"""
    close_logger(logger::Logger)

Close all open file handles in the logger.
"""
function close_logger(logger::Logger)
    if logger.file_handle !== nothing
        close(logger.file_handle)
    end
    
    if logger.state_file_handle !== nothing
        close(logger.state_file_handle)
    end
end

"""
    log_message(logger::Logger, level::Symbol, message::String)

Log a message with the specified level.
"""
function log_message(logger::Logger, level::Symbol, message::String)
    # Define level priorities
    level_priority = Dict(
        :debug => 1,
        :info => 2,
        :warn => 3,
        :error => 4
    )
    
    # Only log messages at or above the current log level
    if get(level_priority, level, 0) >= get(level_priority, logger.log_level, 0)
        timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
        log_line = "[$timestamp][$level] $message"
        
        # Write to file if available
        if logger.file_handle !== nothing
            println(logger.file_handle, log_line)
            flush(logger.file_handle)
        end
        
        # Write to console if enabled
        if logger.console_output
            println(log_line)
        end
    end
end

"""
    log_info(logger::Logger, message::String)

Log an info-level message.
"""
function log_info(logger::Logger, message::String)
    log_message(logger, :info, message)
end

"""
    log_warn(logger::Logger, message::String)

Log a warning-level message.
"""
function log_warn(logger::Logger, message::String)
    log_message(logger, :warn, message)
end

"""
    log_error(logger::Logger, message::String)

Log an error-level message.
"""
function log_error(logger::Logger, message::String)
    log_message(logger, :error, message)
end

"""
    log_event(logger::Logger, time::Float64, event::String)

Log a simulation event with timestamp.
"""
function log_event(logger::Logger, time::Float64, event::String)
    log_info(logger, "t=$(round(time, digits=3)): $event")
end

"""
    log_states(logger::Logger, time::Float64, states::AbstractVector{Float64}, indices::StateIndices)

Log selected state variables at the given simulation time.
"""
function log_states(logger::Logger, time::Float64, states::AbstractVector{Float64}, indices::StateIndices)
    # Only log if we have a state file handle
    if logger.state_file_handle === nothing
        return
    end
    
    # Round time to avoid floating point issues
    time_str = string(round(time, digits=6))
    
    # Log machine states
    log_state_value(logger, time_str, "delta", states[indices.machine[:DELTA]])
    log_state_value(logger, time_str, "omega", states[indices.machine[:OMEGA]])
    
    # Log network states - bus voltages
    log_state_value(logger, time_str, "v1_d", states[indices.network[:V_1_D_IDX]])
    log_state_value(logger, time_str, "v1_q", states[indices.network[:V_1_Q_IDX]])
    log_state_value(logger, time_str, "v2_d", states[indices.network[:V_2_D_IDX]])
    log_state_value(logger, time_str, "v2_q", states[indices.network[:V_2_Q_IDX]])
    log_state_value(logger, time_str, "v3_d", states[indices.network[:V_3_D_IDX]])
    log_state_value(logger, time_str, "v3_q", states[indices.network[:V_3_Q_IDX]])
    
    # Log AVR states
    log_state_value(logger, time_str, "efd", states[indices.avr[:EFD_IDX]])
    
    # Log governor states
    log_state_value(logger, time_str, "fv", states[indices.governor[:FV_IDX]])
end

"""
    log_state_value(logger::Logger, time::String, var_name::String, value::Float64)

Write a single state variable value to the state log file.
"""
function log_state_value(logger::Logger, time::String, var_name::String, value::Float64)
    println(logger.state_file_handle, "$time,$var_name,$value")
end

end # module