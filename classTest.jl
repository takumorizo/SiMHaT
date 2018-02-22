module __myClass
    # export すると, using myClass; myClass.testPrivate() で関数を call 出来てしまう.
    # この状態であれば, myClass.testPrivate(), testPrivate() は call 出来ない
    function testPrivate()::Void
        println("this is a private function!!")
    end
end

module myClass
    export MyClass
    using __myClass

    export MyClass
    type MyClass{T}
        x::T
    end

    function init(x::T)::MyClass{T} where T <:Integer
        return MyClass(x)
    end

    function add!(x::MyClass{T})::MyClass{T} where T <:Integer
        __myClass.testPrivate()
        x.x = x.x+1
        return x
    end
end
