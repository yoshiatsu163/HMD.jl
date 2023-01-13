using Test

function test()

#@testset "HierarchyLabels" begin
    l1, l2, l3 = Label(1, "l1"), Label(2, "l2"), Label(3, "l3")
    noexist    = Label(-1, "not_exist_in_h")
    handmade = begin
        lh = LabelHierarchy()
        @assert add_vertex!(_hierarchy(lh)); set_prop!(_hierarchy(lh), 1, :label, l1)
        @assert add_vertex!(_hierarchy(lh)); set_prop!(_hierarchy(lh), 2, :label, l2)
        @assert add_vertex!(_hierarchy(lh)); set_prop!(_hierarchy(lh), 3, :label, l3)
        add_edge!(_hierarchy(lh), 2, 1)
        add_edge!(_hierarchy(lh), 3, 2)
        _label2node(lh)[l1] = 1
        _label2node(lh)[l2] = 2
        _label2node(lh)[l3] = 3
        lh
    end

    # addの順番違いとinsert，エラー発生時の不変性をテスト
    addrel = begin
        lh = LabelHierarchy()
        _add_label!(lh, l1; super = No_Label, sub = No_Label, insert=false)
        _add_label!(lh, l2; super = No_Label, sub = No_Label)
        _add_label!(lh, l3; super = No_Label, sub = No_Label)
        _add_relation!(lh; super = l1, sub = l2)
        _add_relation!(lh; super = l2, sub = l3)
        lh
    end
    addlb1 = begin
        lh = LabelHierarchy()
        _add_label!(lh, l1; super = No_Label, sub = No_Label)
        _add_label!(lh, l2; super = l1      , sub = No_Label)
        _add_label!(lh, l3; super = l2      , sub = No_Label)
        lh
    end
    addlb2 = begin
        lh = LabelHierarchy()
        _add_label!(lh, l1; super = No_Label, sub = No_Label)
        _add_label!(lh, l3; super = No_Label, sub = No_Label)
        _add_label!(lh, l2; super = l1      , sub = l3      )
        lh
    end
    inserted = begin
        lh = LabelHierarchy()
        _add_label!(lh, l1; super = No_Label, sub = No_Label)
        _add_label!(lh, l3; super = l1      , sub = No_Label)
        lh_copy = deepcopy(lh)
        @test_throws ErrorException _add_label!(lh, l2; super = l1      , sub = l3)
        @test lh == lh_copy
        _add_label!(lh, l2; super = l1      , sub = l3      , insert=true)
        lh
    end
    @test handmade == addrel == addlb1 == addlb2 == inserted

    #@test _hierarchy(addlb) == _hierarchy(addrel) == _hierarchy(handmade)
    #@test _label2node(addlb) == _label2node(addrel) == _label2node(handmade)
    #@test _contains(handmade, l1)
    #@test !_contains(handmade, noexist)
    #@test _get_nodeid(handmade, l1) == 1 && _get_nodeid(h, l2) == 2 && _get_nodeid(h, l3) == 3
    #@test_throws KeyError _get_nodeid(handmade, noexist)
    #@test _super(addlb, l3) == [l2] && _super(addlb, l2) == [l1]
    #@test _sub(addlb, l1)   == [l2] && _sub(addlb, l2)   == [l3]
    #@test _issuper(handmade, l1, l2) && _issuper(handmade, l2, l3)
    #@test _issub(handmade, l2, l1)   && _issub(handmade, l3, l2)

    #@test_throws ErrorException _add_label!(handmade, l1, super = No_Label, sub = No_Label)
    #@test_throws ErrorException _add_label!(handmade, Label(43256, "some"), super = noexist, sub = No_Label)
    #@test_throws ErrorException _add_label!(handmade, Label(43256, "some"), super = No_Label, sub = noexist)
    #@test_throws ErrorException _add_label!(handmade, Label(43256, "some"), super = noexist, sub = noexist)
    #@test_throws ErrorException _add_label!(handmade, Label(43256, "some"), super = l1, sub = l1)
    #@test _hierarchy(handmade)   == _hierarchy(addlb)
    #@test _label2node(handmade) == _label2node(addlb)

    #@test_throws ErrorException _add_relation!(handmade, super = noexist, sub = l1)
    #@test_throws ErrorException _add_relation!(handmade, super = l1, sub = noexist)
    #@test_throws ErrorException _add_relation!(handmade, super = l1, sub = l1)
    #@test_throws ErrorException _add_relation!(handmade, super = l1, sub = l2)
    #@test_throws ErrorException _add_relation!(handmade, super = l2, sub = l1)
    #@test _hierarchy(handmade)   == _hierarchy(addlb)
    #@test _label2node(handmade) == _label2node(addlb)





#end

end #test
