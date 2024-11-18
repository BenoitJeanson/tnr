function create_case(case::String)
    sys = System(joinpath("data","exp_raw","case", "$case.m"));
    coord = load_coord(joinpath("data","exp_raw","coord", "$case.csv"))
    g = network2graph(sys)
    balance!(g)
    add_constraint(g, b->b.p_max=1)
    bus_confs = [BusConf(6, [SubBus(-.17 , [7, 9])])]
    g, bus_confs, coord
end

function store_result(model, filename)
    open(filename, "w") do file
        foreach(k->println(file, "$k, $(value(k))"), all_variables(model))
    end
end

