include("Include.jl")

# define some colors -
const shaded_color_value = "lightgray"
const mean_color_value = "dimgray"
const experimental_color_value = "black"

const P1_color = "blue"
const P1_shaded_color="skyblue"

const P2_color = "orange"
const P2_shaded_color="navajowhite"

const P3_color = "green"
const P3_shaded_color="lightgreen"

# Script to solve the balance equations -
time_start = 0.0
time_stop = 240.0
time_step_size = 0.1
number_of_timesteps = length(time_start:time_step_size:time_stop)

# How many samples do we want to explore?
number_of_samples = 10
sigma = 0.20

# initialize storage -
time_array = []
P1_array = zeros(number_of_timesteps,number_of_samples)
P2_array = zeros(number_of_timesteps,number_of_samples)
P3_array = zeros(number_of_timesteps,number_of_samples)

# main loop -
for sample_index = 1:number_of_samples

  # Load the data dictionary (default parameter values) -
  data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

  # no connection between P3 and P2 -
  control_parameter_dictionary = data_dictionary["control_parameter_dictionary"];
  control_parameter_dictionary["W_gene_2_gene_3"] = 0.0;
  data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary


  # Pertub the RNAPs and Ribosome abundance -
  rnapII_concentration  = data_dictionary["rnapII_concentration"] # muM
	ribosome_concentration  = data_dictionary["ribosome_concentration"] # muM

  rnapII_concentration = abs(rnapII_concentration*(1+sigma*randn()))
  ribosome_concentration  = abs(ribosome_concentration*(1+sigma*randn()))

  data_dictionary["rnapII_concentration"] = rnapII_concentration  # muM
	data_dictionary["ribosome_concentration"] = ribosome_concentration # muM

  # Run the washout simulation -
  (T,X) = washout_simulation(time_start,time_stop,time_step_size,data_dictionary)

  time_array = T
  for (time_index,time_value) in enumerate(T)
    P1_array[time_index,sample_index] = X[time_index,7];
    P2_array[time_index,sample_index] = X[time_index,8];
    P3_array[time_index,sample_index] = X[time_index,9];
  end
end

# Confidence interval?
SF = (1.96/sqrt(number_of_samples))

# Make some plots -
# P1 -
P1_mean = mean(P1_array,2)
P1_std = std(P1_array,2)
P1_lower_bound = P1_mean - SF*P1_std
P1_upper_bound = P1_mean + SF*P1_std
plot(time_array,P1_mean,lw=2,color=P1_color)
fill_between(time_array,vec(P1_lower_bound),vec(P1_upper_bound),color=P1_shaded_color,lw=3)

# P2 -
P2_mean = mean(P2_array,2)
P2_std = std(P2_array,2)
P2_lower_bound = P2_mean - SF*P2_std
P2_upper_bound = P2_mean + SF*P2_std
plot(time_array,P2_mean,lw=2,color=P2_color)
fill_between(time_array,vec(P2_lower_bound),vec(P2_upper_bound),color=P2_shaded_color,lw=3)

# P3 -
P3_mean = mean(P3_array,2)
P3_std = std(P3_array,2)
P3_lower_bound = P3_mean - SF*P3_std
P3_upper_bound = P3_mean + SF*P3_std
plot(time_array,P3_mean,lw=2,color=P3_color)
fill_between(time_array,vec(P3_lower_bound),vec(P3_upper_bound),color=P3_shaded_color,lw=3)

# Dump to disk -
savefig("../figs/Washout-3G-KO.pdf")
