using Test

function test()
    
@testset "heierchy" begin
    h = Heierchy()
    l1 = Label(1, "l1")
    @test add_vertex!(_heierchy(h))
    set_prop!(_heierchy(h), 1, :label, l1)

    @test _get_nodeid(h, l1; allow_NA = false) == 1
    @test _get_nodeid(h, l1; allow_NA = true)  == 1
    @test _get_nodeid(h, Label(2, "l1"); allow_NA = true) |> ismepty
    @test _get_nodeid(h, Label(1, "l2"); allow_NA = true) |> ismepty
    @test_throws ErrorException _get_nodeid(h, Label(2, "l1"); allow_NA = false)
    @test_throws ErrorException _get_nodeid(h, Label(1, "l2"); allow_NA = false)

    hh = Heierchy()
    _add_label!(hh, l1; super = nothing, sub = nothing)
    @test h == hh
    @test_throws ErrorException _add_label!(hh, l1; super = Label(2, "l2"), sub = nothing)
    @test_throws ErrorException _add_label!(hh, l1; super = nothing, sub = Label(2, "l2"))
    
    l2 = Label(2, "l2")
    _add_label!(hh, l2, super = l1, sub = nothing)
    @test_throws ErrorException _add_label!(hh, l1; super = l1, sub = l2)
    current_id = _get_nodeid(hh, l2; allow_NA = false)
    super_id   = _get_nodeid(hh, l1; allow_NA = false)
    @test current_id == 2
    @test super_id   == 1
    @test _super(hh, l2) == [super_id]
    @test _sub(hh, l1)   == [current_id]
    @test _issuper(hh, l1, l2)
    @test _issub(hh, l2, l1)

    #add_relation
end

end #test