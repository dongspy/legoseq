use minijinja::{
    context,
    value::{Object, Value},
    Environment, Template,
};
use serde::{Deserialize, Serialize};
use std::{
    collections::HashMap,
    fmt::{self, Display},
    fs,
    hash::Hash,
};

#[derive(Debug, Deserialize, Serialize)]
struct SeqOut {
    name: String,
    seq: Option<String>,
    qual: Option<String>,
}
impl Display for SeqOut {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Use `self.number` to refer to each positional data point.
        write!(f, "{}", self.name)
    }
}

impl Object for SeqOut {}

#[test]
fn test() {
    // 读取模板文件
    let template_string = fs::read_to_string("./test/template.txt").expect("无法读取模板文件");
    let mut seq_hash = HashMap::new();
    let seqout1 = SeqOut {
        name: "name".to_string(),
        seq: Some("seq".to_string()),
        qual: None,
    };
    seq_hash.insert("aa", seqout1);
    let ctx = Value::from_serializable(&seq_hash);
    // 创建一个新的 MiniJinja 环境
    let env = Environment::new();

    // 从字符串创建一个模板
    let template = env
        .template_from_str(&template_string)
        .expect("无法从字符串创建模板");

    // let ctx: minijinja::value::Value = context!(variable_name);
    // 渲染模板
    let result = template.render(ctx).expect("无法渲染模板");

    println!("{}", result);
}
