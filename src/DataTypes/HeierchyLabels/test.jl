using Test

function test()

@testset "HeierchyLabels" begin
    l1, l2, l3 = Label(1, "l1"), Label(2, "l2"), Label(3, "l3")
    noexist    = Label(-1, "not_exist_in_h")
    @testset "heierchy" begin
        handmade = begin
            h = Heierchy()
            @assert add_vertex!(_heierchy(h)); set_prop!(_heierchy(h), 1, :label, l1)
            @assert add_vertex!(_heierchy(h)); set_prop!(_heierchy(h), 2, :label, l2)
            @assert add_vertex!(_heierchy(h)); set_prop!(_heierchy(h), 3, :label, l3)
            add_edge!(_heierchy(h), 2, 1)
            add_edge!(_heierchy(h), 3, 2)
            _label2node(h)[l1] = 1
            _label2node(h)[l2] = 2
            _label2node(h)[l3] = 3
            h
        end
        addlb = begin
            h = Heierchy()
            _add_label!(h, l1; super = No_Label, sub = No_Label)
            _add_label!(h, l2; super = l1     , sub = No_Label)
            _add_label!(h, l3; super = l2     , sub = No_Label)
            h
        end
        addrel = begin
            h = Heierchy()
            _add_label!(h, l1; super = No_Label, sub = No_Label)
            _add_label!(h, l2; super = No_Label, sub = No_Label)
            _add_label!(h, l3; super = No_Label, sub = No_Label)
            _add_relation!(h, super = l1, sub = l2)
            _add_relation!(h, super = l2, sub = l3)
            h
        end

        @test _heierchy(addlb) == _heierchy(addrel) == _heierchy(handmade)
        @test _label2node(addlb) == _label2node(addrel) == _label2node(handmade)
        @test _contains(handmade, l1)
        @test !_contains(handmade, noexist)
        @test _get_nodeid(handmade, l1) == 1 && _get_nodeid(h, l2) == 2 && _get_nodeid(h, l3) == 3
        @test_throws KeyError _get_nodeid(handmade, noexist)
        @test _super(addlb, l3) == [l2] && _super(addlb, l2) == [l1]
        @test _sub(addlb, l1)   == [l2] && _sub(addlb, l2)   == [l3]
        @test _issuper(handmade, l1, l2) && _issuper(handmade, l2, l3)
        @test _issub(handmade, l2, l1)   && _issub(handmade, l3, l2)

        @test_throws ErrorException _add_label!(handmade, l1, super = No_Label, sub = No_Label)
        @test_throws ErrorException _add_label!(handmade, Label(43256, "some"), super = noexist, sub = No_Label)
        @test_throws ErrorException _add_label!(handmade, Label(43256, "some"), super = No_Label, sub = noexist)
        @test_throws ErrorException _add_label!(handmade, Label(43256, "some"), super = noexist, sub = noexist)
        @test_throws ErrorException _add_label!(handmade, Label(43256, "some"), super = l1, sub = l1)
        @test _heierchy(handmade)   == _heierchy(addlb)
        @test _label2node(handmade) == _label2node(addlb)

        @test_throws ErrorException _add_relation!(handmade, super = noexist, sub = l1)
        @test_throws ErrorException _add_relation!(handmade, super = l1, sub = noexist)
        @test_throws ErrorException _add_relation!(handmade, super = l1, sub = l1)
        @test_throws ErrorException _add_relation!(handmade, super = l1, sub = l2)
        @test_throws ErrorException _add_relation!(handmade, super = l2, sub = l1)
        @test _heierchy(handmade)   == _heierchy(addlb)
        @test _label2node(handmade) == _label2node(addlb)
    end

    @testset "label2atom" begin
        handmade = begin
            l2a = Label2atom()
            push!(l2a.toatom, l1 => Set([1,2,3,4]))
            for i in [1,2,3,4]
                push!(l2a.tolabel, i => Set([l1]))
            end
            l2a
        end
        addlb = begin
            l2a = Label2atom()
            _add_atoms!(l2a, 4)
            _add_label!(l2a, l1, [1,2,3,4])
            l2a
        end

        @test tolabel(handmade) == tolabel(addlb)
        @test toatom(handmade)  == toatom(addlb)
        @test _labels(addlb) == Set([l1])
        @test addlb[l1] == Set([1,2,3,4])
        @test addlb[1] == addlb[2] == addlb[3] == addlb[4] == Set([l1])

        @test_throws ErrorException _add_label!(addlb, l1, [1,2,3,4])
        @test_throws KeyError addlb[l2]
        @test_throws KeyError addlb[6]
        @test_throws ErrorException begin
            tmp = Label2atom()
            _add_label!(addlb, l1, [1,2,3,4])
        end
        @test tolabel(handmade) == tolabel(addlb)
        @test toatom(handmade)  == toatom(addlb)
    end

    addrel = begin
        addrel = LabelHeierchy()
        add_atoms!(addrel, 4)
        add_label!(addrel, l1, [1,2,3,4]; super = No_Label, sub = No_Label)
        add_label!(addrel, l2, [1,2,3]  ; super = No_Label, sub = No_Label)
        add_label!(addrel, l3, [1,2]    ; super = No_Label, sub = No_Label)
        add_relation!(addrel; super = l1, sub = l2)
        add_relation!(addrel; super = l2, sub = l3)
        addrel
    end

    addlb1 = begin
        addlb = LabelHeierchy()
        add_atoms!(addlb, 4)
        add_label!(addlb, l1, [1,2,3,4]; super = No_Label, sub = No_Label)
        add_label!(addlb, l2, [1,2,3]  ; super = l1     , sub = No_Label)
        add_label!(addlb, l3, [1,2]    ; super = l2     , sub = No_Label)
        addlb
    end
    addlb2 = begin
        addlb = LabelHeierchy()
        add_atoms!(addlb, 4)
        add_label!(addlb, l1, [1,2,3,4]; super = No_Label, sub = No_Label)
        add_label!(addlb, l3, [1,2]    ; super = No_Label, sub = No_Label)
        add_label!(addlb, l2, [1,2,3]  ; super = l1     , sub = l3)
        addlb
    end

    #一致確認
    #@test addlb1 == addlb2 == addrel
    @test addlb1 == addrel
end

end #test
