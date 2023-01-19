using Test

function test()

@testset "HierarchyLabels" begin
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
    addrel1 = begin
        lh = LabelHierarchy()
        _add_label!(lh, l1)
        _add_label!(lh, l2)
        _add_label!(lh, l3)
        _add_relation!(lh; super = l1, sub = l2)
        _add_relation!(lh; super = l2, sub = l3)
        lh
    end
    addrel2 = begin
        lh = LabelHierarchy()
        _add_label!(lh, l1)
        _add_label!(lh, l2)
        _add_label!(lh, l3)
        _add_relation!(lh; super = l2, sub = l3)
        _add_relation!(lh; super = l1, sub = l2)
        lh
    end

    @test _hierarchy(addrel1) == _hierarchy(addrel2) == _hierarchy(handmade)
    @test _contains(addrel1, l1)
    @test !_contains(addrel1, noexist)
    @test _get_nodeid(addrel1, l1) == 1 && _get_nodeid(addrel1, l2) == 2 && _get_nodeid(addrel1, l3) == 3
    @test _super(addrel1, l3) == [l2] && _super(addrel1, l2) == [l1]
    @test _sub(addrel1, l1)   == [l2] && _sub(addrel1, l2)   == [l3]
    @test _issuper(addrel1, l1, l2) && _issuper(addrel1, l2, l3)
    @test _issub(addrel1, l2, l1)   && _issub(addrel1, l3, l2)

    @test_throws KeyError _get_nodeid(addrel1, noexist)
    @test_throws ErrorException _has_relation(addrel1, noexist, l2)
    @test_throws ErrorException _get_label(LabelHierarchy(), 0)
    @test_throws KeyError _get_label(addrel1, -1)

    @test _add_label!(addrel1, l1) == Label_Occupied
    @test _add_relation!(addrel1; super=noexist, sub=l3) == Label_Missing
    @test _add_relation!(addrel1; super=l3, sub=l3)      == Label_Duplication
    @test _add_relation!(addrel1; super=l2, sub=l3)      == Relation_Occupied

    lh = deepcopy(addrel1)
    @test lh == addrel1
    @test_throws KeyError _remove_relation!(lh, noexist, l3)
    @test_throws KeyError _remove_relation!(lh, l2, noexist)
    @test lh == addrel1
    @test _remove_label!(lh, l1)
    @test _remove_relation!(lh, l2, l3)
    _remove_label!(lh, l1)
end

end #test
