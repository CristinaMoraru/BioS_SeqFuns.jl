export terminal_repeat

@kwdef struct terminal_repeat
    type::String
    length::Int
    left_start::Int = 1
    left_end::Int
    right_start::Int
    right_end::Int
end