

def switch_class(classes,to_add,all_classes):
    classes = set(classes.split(" "))
    for c in all_classes:
        if c in classes:
            classes.remove(c)
    for c in to_add:
        classes.add(c)
    
    return " ".join(classes)